#' Liftover a GRanges from one genome build to another
#'
#' @param genomic_region A GRanges object.
#' @param chain A "Chain" object to be used with rtracklayer::liftOver.
#' @param remove_non_mapping A logical value indicating whether to remove regions in the source genome.
#' which could not be mapped to the target genome. Default is TRUE. 
#' @param remove_one_to_many_mapping A logical value indicating whether to remove regions in the source genome.
#' which map to multiple regions in the target genome. Default is TRUE.
#' @param remove_many_to_one_mapping A logical value indicating whether to remove different regions in the source genome.
#' which map to the same region in the target genome i.e. remove regions which overlap after the liftover. Default is TRUE.
#' @param target_regions_overlap_filter An optional GRanges object used to filter the rowRanges by overlaps after liftover, 
#' for example CpG sites from the target genome. Any regions which do not overlap target_regions_overlap_filter will be removed.  
#' @param flatten_ranges A logical value indicating whether to convert the lifted over ranges from a 
#' GRangesList to GRanges if all remaining source regions can be uniquely mapped to the target genome. 
#' @return A GRanges object
#' @export
liftover_granges = function(genomic_regions, chain, remove_non_mapping = T, remove_one_to_many_mapping = T, 
  remove_many_to_one_mapping = T, target_regions_overlap_filter = NULL, flatten_ranges = T){
  
  # Liftover genomic_regions using specified liftover chain file
  liftover_ranges = rtracklayer::liftOver(genomic_regions, chain)
  
  # Put seqlevels of liftover_ranges in the same order as genomic_region
  seqlevels(liftover_ranges) = seqlevels(genomic_regions)
  
  # Initialize selected regions to all liftover_ranges
  selected_ranges = seq_along(liftover_ranges)
  
  # Get the number of regions in the target genome each region in the source genome matches to
  mappings_count = lengths(liftover_ranges)
  
  # Remove non-mapping regions from selected_ranges if specified
  if(remove_non_mapping){
    non_mapping_regions = which(mappings_count == 0)
    selected_ranges = setdiff(selected_ranges, non_mapping_regions)
  }
  
  # Remove one-to-many mapping regions from selected_ranges if specified
  if(remove_one_to_many_mapping){
    one_to_many_mapping_regions = which(mappings_count > 1)
    selected_ranges = setdiff(selected_ranges, one_to_many_mapping_regions)
  }
  
  # Remove many-to-one mapping regions from selected_ranges if specified
  if(remove_many_to_one_mapping){
    self_overlaps = countOverlaps(liftover_ranges, liftover_ranges)
    many_to_one_mapping_regions = which(self_overlaps > 1)
    selected_ranges = setdiff(selected_ranges, many_to_one_mapping_regions)
  }
  
  # Identify regions which overlap target_regions_overlap_filter is provided
  if(!is.null(target_regions_overlap_filter)){
    target_regions_overlaps = countOverlaps(liftover_ranges, target_regions_overlap_filter)
    non_target_overlaps = which(target_regions_overlaps < 1)
    selected_ranges = setdiff(selected_ranges, non_target_overlaps)
  }
  
  # Subset liftover_ranges for selected ranges
  liftover_ranges = liftover_ranges[selected_ranges]
  
  # Flatten liftover_ranges if specified and there are no one-to-many mappings
  if(flatten_ranges & length(unlist(liftover_ranges)) == length(liftover_ranges)){
    liftover_ranges = unlist(liftover_ranges)
  }
  
  # Return liftover_ranges 
  return(liftover_ranges)
  
}

#' Liftover bigWig or bedGraph files from one genome build to another  
#'
#' Liftover either bigwig or bedGraphs from hg19 to hg38 or from hg38 to hg19 and save as a bigWig or a bedGraph. 
#' bedRemoveOverlap is used to remove overlapping regions resulting from the output files after liftover. 
#' Output files have the basename of the input files with "liftover" and either "hg19" or "hg38" appended, depending on value of to_genome_build 
#'
#' @param input_filepaths Paths to either bedGraph files or bigWig files to be lifted-over. All files must be of the same type. bedGraph files may be zipped. 
#' @param input_filetype Either "bw" for bigWig or "bg" for bedGraph.
#' @param output_filetype Either "bw" for bigWig or "bg" for bedGraph.
#' @param chain A "Chain" object to be used with rtracklayer::liftOver.
#' @param to_genome_build Either "hg19" or "hg38" depending on the intended output genome
#' @param output_directory Directory to save output files in. Default is current working directory. If output_directory doesn't exist, it is created.
#' @param compress_output_bedGraphs A logical value indicating whether output bedGraph files should be compressed. Default is TRUE. 
#' @param ncores The number of cores to use to process files in parallel. Default is 1.
#' @return None
#' @export
liftover_files = function(input_filepaths, input_filetype, output_filetype, to_genome_build, output_directory = ".", compress_output_bedGraphs = T, ncores = 1){
  
  # Check that allowed values are entered for input_filetype, output_filetype, to_genome_build and ncores
  match.arg(arg = input_filetype, choices = c("bw", "bg"), several.ok = F)
  match.arg(arg = output_filetype, choices = c("bw", "bg"), several.ok = F)
  match.arg(arg = to_genome_build, choices = c("hg19", "hg38"), several.ok = F)
  if(ncores %% 1 != 0 | ncores < 0){stop("ncores must be a positive whole number")}
  
  # Create output_directory if it doesn't exist)
  dir.create(output_directory, recursive = T, showWarnings = F)
  
  # Make a cluster 
  if(ncores > 1){
    cluster = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cluster, cores = ncores)
  }
  `%dopar%` = foreach::`%dopar%`
  
  # Create output filenames 
  output_filenames = paste(output_directory, paste0(tools::file_path_sans_ext(gsub(".gz", "", input_filepaths)), "_liftover_", to_genome_build, ".", tools::file_ext(input_filepaths)), sep = "/")
  
  # Load correct chain file for relevent genome build
  if (to_genome_build == "hg19"){
    chain_file = system.file("extdata", "hg38ToHg19.over.chain", package = "genomicTools")
    chrom_sizes = system.file("extdata", "hg19_chrom_size.tsv", package = "genomicTools")
  } else if (to_genome_build == "hg38"){
    chain_file = system.file("extdata", "hg19ToHg38.over.chain", package = "genomicTools")
    chrom_sizes = system.file("extdata", "hg38_chrom_size.tsv", package = "genomicTools")
  }
  
  # Liftover files
  foreach::foreach(file_number = seq_along(input_filepaths)) %dopar% {
  
    # Set current working file to the input file
    current_file = input_filepaths[file_number]
    
    # If input type is bw, convert to a bedGraph
    system("echo Converting to bedGraph")
    bg_temp = tempfile()
    if(input_filetype == "bw"){
      system(paste(system.file("extdata", "bigWigToBedGraph", package = "genomicTools"), input_filepaths[file_number], bg_temp))
      current_file = bg_temp
    }
    
    # LiftOver bedGraph
    liftover_temp = tempfile()
    system(paste(system.file("extdata", "liftOver", package = "genomicTools"), current_file, chain_file, liftover_temp, "/dev/null"))
    current_file = liftover_temp
    unlink(bg_temp)
    
    # Remove duplicated regions post liftover and sort
    system("echo Removing duplicated regions")
    overlap_removed_temp = tempfile()
    system(paste("sort -k1,1 -k2,2n -o", current_file, current_file))
    system(paste(system.file("extdata", "bedRemoveOverlap", package = "genomicTools"), current_file, overlap_removed_temp))
    current_file = overlap_removed_temp
    unlink(liftover_temp)
    
    # If output_filetype is bg, can return current_file. If output_filetype is bw, need to convert it back to a bigwig
    if(output_filetype == "bg"){
      system(paste("mv", current_file, output_filenames))
    } else if(output_filetype == "bw"){
      system("echo Converting to bigWig")
      system(paste(system.file("extdata", "bedGraphToBigWig", package = "genomicTools"), current_file, chrom_sizes, output_filenames[file_number]))
    }
    
    # Gzip output bedGraph if specified
    if(output_filetype == "bg"){
      if(compress_output_bedGraphs){
      system(paste("gzip", output_filenames))
      }
    }
    
    # Remove temporary files
    unlink(overlap_removed_temp)
  }
  
}