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
    chain_file = system.file("extdata", "hg38ToHg19.over.chain", package = "genomeTools")
    chrom_sizes = system.file("extdata", "hg19_chrom_size.tsv", package = "genomeTools")
  } else if (to_genome_build == "hg38"){
    chain_file = system.file("extdata", "hg19ToHg38.over.chain", package = "genomeTools")
    chrom_sizes = system.file("extdata", "hg38_chrom_size.tsv", package = "genomeTools")
  }
  
  # Liftover files
  foreach::foreach(file_number = seq_along(input_filepaths)) %dopar% {
  
    # Set current working file to the input file
    current_file = input_filepaths[file_number]
    
    # If input type is bw, convert to a bedGraph
    system("echo Converting to bedGraph")
    bg_temp = tempfile()
    if(input_filetype == "bw"){
      system(paste(system.file("extdata", "bigWigToBedGraph", package = "genomeTools"), input_filepaths[file_number], bg_temp))
      current_file = bg_temp
    }
    
    # LiftOver bedGraph
    liftover_temp = tempfile()
    system(paste(system.file("extdata", "liftOver", package = "genomeTools"), current_file, chain_file, liftover_temp, "/dev/null"))
    current_file = liftover_temp
    unlink(bg_temp)
    
    # Remove duplicated regions post liftover and sort
    system("echo Removing duplicated regions")
    overlap_removed_temp = tempfile()
    system(paste("sort -k1,1 -k2,2n -o", current_file, current_file))
    system(paste(system.file("extdata", "bedRemoveOverlap", package = "genomeTools"), current_file, overlap_removed_temp))
    current_file = overlap_removed_temp
    unlink(liftover_temp)
    
    # If output_filetype is bg, can return current_file. If output_filetype is bw, need to convert it back to a bigwig
    if(output_filetype == "bg"){
      system(paste("mv", current_file, output_filenames))
    } else if(output_filetype == "bw"){
      system("echo Converting to bigWig")
      system(paste(system.file("extdata", "bedGraphToBigWig", package = "genomeTools"), current_file, chrom_sizes, output_filenames[file_number]))
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