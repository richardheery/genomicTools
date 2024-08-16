#' Summarize values in bigWig files for a set of genomic regions
#' 
#' Calculates the mean (with or without non-covered bases counting as zeros) or the sum of values in bigWig files 
#' for provided regions. Regions can be provided either as a GRanges object or as a path to a BED file.
#' Columns in the output correspond to input bigWig files and rows correspond to genomic regions. 
#' Uses names of the provided BED file or GRanges as row names so long as they are unique. 
#' Otherwise names regions region_1, region_2, etc. 
#'
#' @param bw_filepaths A vector of filepaths for bigWig files.
#' @param bed_filepath Filepath for a BED file with genomic regions of interest. 
#' One of either bed_filepath or grs must be provided. 
#' @param gr A GRanges object. One of either bed_filepath or gr must be provided. 
#' @param statistic Statistic to compute. One of "mean" (average in bigWig files over just covered bases in regions), 
#' "mean0" (average over bases with non-covered bases counting as zeroes) or "sum".
#' @param column_names A vector of names to use as column names of the output table. 
#' Default is to use the base name of bw_filepaths.
#' @param ncores Number of cores to use. Default is 1
#' @return A data.frame with the results.
#' @export
bigwig_summarize_over_regions = function(bw_filepaths, bed_filepath = NULL, 
  gr = NULL, statistic, column_names = NULL, ncores = 1){
  
  # Check that one of bed_filepath or gr is provided
  if(!is.null(gr)){
    if(!is.null(bed_filepath)){
      stop("Either bed_filepath or gr should be provided, not both")
    } else {
      # If gr is provided, remove the metadata, and export it as a temporary BED file 
      mcols(gr) = NULL
      bed_filepath = tempfile(pattern = "temp_bed")
      rtracklayer::export.bed(gr, bed_filepath)
    }
  } else if(is.null(bed_filepath)){
    # Throw an error if neither bed_filepath or gr are provided
    stop("One of bed_filepath or gr must be provided")
  }
  
  # Check that one of the correct choices is provided for statistic
  match.arg(arg = statistic, choices = c("mean", "mean0", "sum"), several.ok = F)
  
  # Check that if column_names is provided, it has the same length as bw_filepaths
  if(!is.null(column_names)){if(length(bw_filepaths) != length(column_names)){stop("Length of bw_filepaths not equal to length of column_names")}}
  
  # Get names of regions from BED file and check that they are unique
  region_names = data.table::fread(bed_filepath, header = F, sep = "\t", stringsAsFactors = F, select = 4)$V4
  if(any(duplicated(region_names))){
    message("Names of regions are not unique. Will thus name them region_1, region_2, etc.")
    region_names = paste("region", seq_along(region_names))
  }
  
  # Select the mean, mean0 or sum  columns from the bigWigAverageOverBed output files for downstream use
  statistic_col = ifelse(statistic == "sum", 4, ifelse(statistic == "mean0", 5, 6))
  
  # Create temporary output directory
  temp_dir = tempfile("bigwig_summarize_over_regions_tmp/")
  dir.create(temp_dir)
  
  # Create names for the output files that will be produced by bigWigAverageOverBed  
  summary_over_bed_file_names = paste0(tools::file_path_sans_ext(basename(bw_filepaths)), "_summary_over_bed.txt")
  
  # If no sample names provided, set to the basename of the input bigwig-files
  if(is.null(column_names)){column_names = basename(tools::file_path_sans_ext(bw_filepaths))}
  
  # Make a cluster if more than one core being used
  if(ncores > 1){
    cluster = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cluster, cores = ncores)
    on.exit(parallel::stopCluster(cluster))
    `%dopar%` = foreach::`%dopar%`
  } else {
    `%dopar%` = foreach::`%do%`
  }
  
  # Use bigWigAverageOverBed to calculate mean value for each bigwig file for each region in the BED file and save the results in output directory
  foreach::foreach(bigwig_number = seq_along(bw_filepaths)) %dopar% {
    system(paste(system.file("extdata", "bigWigAverageOverBed", package = "genomeTools"), bw_filepaths[bigwig_number], bed_filepath, 
      paste0(temp_dir, summary_over_bed_file_names[bigwig_number])))}
  
  # Create a matrix with rows corresponding to genomic regions and columns corresponding to genomic features
  region_values = data.frame(foreach::foreach(result_file = 
   list.files(temp_dir, full.names = T), .combine = "cbind") %dopar% {
     data.table::fread(result_file, sep = "\t", header = F, select = statistic_col, nThread = 1)
  })
  
  # Remove temporary files
  unlink(temp_dir, recursive = T)
  if(!is.null(gr)){
    unlink(bed_filepath)
  }
  
  # Assign rownames and column names to dataframe
  rownames(region_values) = region_names
  colnames(region_values) = column_names
  region_values = tibble::rownames_to_column(region_values, "genomic_region")
  
  # Return the 
  return(region_values)
  
}
