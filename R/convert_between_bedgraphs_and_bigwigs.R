#' Convert begGraph files to bigWig files 
#'
#' Convert begGraph files to bigWig files and save in a specified directory. 
#' bedGraph files can be compressed with gzip. Takes ~ 1 minute per file. If more than one file, each file can be processed in parallel.
#'
#' @param bg_files Paths to a directory of bedGraph files.
#' @param chrom_sizes_file Path to a file with the chromosome sizes for the relevant genome.
#' @param output_dir Name of the directory where created bigWigs will be saved. Default is current working directory.
#' @param ncores The number of cores to use. Default 1.
#' @return None
#' @export
convert_bedgraphs_to_bigwigs = function(bg_files, chrom_sizes_file, output_dir = ".", ncores = 1){
  
  # Check thatncores is a positive whole integer
  if(ncores %% 1 != 0 | ncores < 0){stop("ncores must be a positive whole number")}
  
  # Create output_dir if it doesn't exist
  dir.create(output_dir, showWarnings = F)
  
  # Make a cluster if more than one core being used
  if(ncores > 1){
    cluster = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cluster, cores = ncores)
  }
  `%do%` = foreach::`%do%`
  `%dopar%` = foreach::`%dopar%`
  
  foreach::foreach(bg = bg_files) %dopar% {
    
    name = basename(gsub(".gz", "", bg))
    
    # If file is compressed, create a temporary uncompressed file
    if(tools::file_ext(bg) == "gz"){
      temp = tempfile(fileext = ".bg")
      system(paste("zcat", bg, ">", temp))
      bg = temp
    }
    
    bg_extension = tools::file_ext(bg)
    output_name = gsub(bg_extension, "bw", name)
    
    system(paste("bedGraphToBigWig", bg, chrom_sizes_file, paste(output_dir, output_name, sep = "/")))
    
    # Remove temporary file
    if(tools::file_ext(bg) == "gz"){unlink(temp)}
    }
    
  }
  