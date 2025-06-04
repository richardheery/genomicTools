#' Convert begGraph files to bigWig files 
#'
#' Convert begGraph files to bigWig files and save in a specified directory. 
#' bedGraph files can be compressed with gzip. Takes ~ 1 minute per file. If more than one file, each file can be processed in parallel.
#'
#' @param bg_files Paths to a directory of bedGraph files 
#' @param genome_build Either "hg19" or "hg38" depending on the intended output genome
#' @param output_dir Name of the directory where created bigWigs will be saved. Default is current working directory
#' @param ncores The number of cores to use. Default 1
#' @return None
#' @export
convert_bedgraphs_to_bigwigs = function(bg_files, genome_build, output_dir = ".", ncores = 1){
  
  # Check that allowed values are entered for genome_build and theat ncores is a positive whole integer
  match.arg(arg = genome_build, choices = c("hg19", "hg38"), several.ok = F)
  if(ncores %% 1 != 0 | ncores < 0){stop("ncores must be a positive whole number")}
  
  # Create output_dir if it doesn't exist
  dir.create(output_dir, showWarnings = F)
  
  # Set correct chromosome sizes file for genome build
  if(genome_build == "hg19"){
    chrom_sizes_file = "~/genomes/genome_files/hg19_chrom_size.tsv"
  } else {
    chrom_sizes_file = "~/genomes/genome_files/hg38_chrom_size.tsv"
  }
  
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
  