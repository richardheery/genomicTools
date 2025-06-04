#' Convert bigWig files to bedGraph files 
#'
#' Convert bigWig files to bedGraph files and save in a specified directory. Can compress bedGraphs using gzip if specified.
#' Takes ~ 1 minute per file. If more than one file, each file can be processed in parallel.
#'
#' @param bw_files Paths to a directory of bigwig files 
#' @param compress_output_bedgraphs A logical value indicating whether to compress output bedGraph files with gzip. Default is TRUE.
#' @param output_dir Name of the directory where created bigWigs will be saved. Default is current working directory
#' @param ncores The number of cores to use. Default 1
#' @return None
#' @export
convert_bigwigs_to_bedgraphs = function(bw_files, compress_output_bedgraphs = T, output_dir = ".", ncores = 1){
  
  # Check that allowed values are entered for genome_build and theat ncores is a positive whole integer
  if(ncores %% 1 != 0 | ncores < 0){stop("ncores must be a positive whole number")}
  
  # Create output_dir if it doesn't exist
  dir.create(output_dir, showWarnings = F)
  
  # Make a cluster if more than one core being used
  if(ncores > 1){
    cluster = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cluster, cores = ncores)
    `%do%` = foreach::`%dopar%`
  } else {
    `%do%` = foreach::`%do%`
  }
  
  foreach::foreach(bw = bw_files) %do% {
    
    # Create output name for bedGraph file by replacing extensionof bigWig file with bg
    output_name = paste0(tools::file_path_sans_ext(basename(bw)), ".bg")
    
    # Convert bigWig to bedGraph file using bigWigToBedGraph and save to output directory
    system(paste(system.file("extdata", "bigWigToBedGraph", package = "genomeTools"), bw,  paste(output_dir, output_name, sep = "/")))
    
    # Compress output directory using gzip if indicated
    if(compress_output_bedgraphs){R.utils::gzip(paste(output_dir, output_name, sep = "/"))}
  }
}
  