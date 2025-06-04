#' #' Liftover either bigwig or bedGraphs 
#' #'
#' #' Liftover either bigwig or bedGraphs from hg19 to hg38 or from hg38 to hg19 and save as a bigWig or a bedGraph
#' #'
#' #' @param input_filepath Path to a bigWig file to be lifted-over 
#' #' @param input_filetype Either "bw" for bigWig or "bg" for bedGraph
#' #' @param output_filetype Either "bw" for bigWig or "bg" for bedGraph
#' #' @param to_genome_build Either "hg19" or "hg38" depending on the intended output genome
#' #' @param output_filename Name to give output file
#' #' @return None
#' #' @export
#' liftover_file_original = function(input_filepath, input_filetype, output_filetype, to_genome_build, output_filename){
#'   
#'   # Check that allowed values are entered for to_genome_build and output_filetype
#'   match.arg(arg = input_filetype, choices = c("bw", "bg"), several.ok = F)
#'   match.arg(arg = output_filetype, choices = c("bw", "bg"), several.ok = F)
#'   match.arg(arg = to_genome_build, choices = c("hg19", "hg38"), several.ok = F)
#'   
#'   # Set current working file to the input file
#'   current_file = input_filepath
#'   
#'   # Load correct chain file for relevent genome build
#'   if (to_genome_build == "hg19"){
#'     chain_file = system.file("extdata", "hg38ToHg19.over.chain", package = "genomeTools")
#'     chrom_sizes = system.file("extdata", "hg19_chrom_size.tsv", package = "genomeTools")
#'   } else if (to_genome_build == "hg38"){
#'     chain_file = system.file("extdata", "hg19ToHg38.over.chain", package = "genomeTools")
#'     chrom_sizes = system.file("extdata", "hg38_chrom_size.tsv", package = "genomeTools")
#'   }
#'   
#'   # Get extension of input file
#'   input_extension = paste0(".", tools::file_ext(input_filepath))
#'   
#'   # If input type is bw, convert to a bedGraph
#'  if(input_filetype == "bw"){
#'     system(paste(system.file("extdata", "bigWigToBedGraph", package = "genomeTools"), input_filepath, basename(gsub(input_extension, "__tempfile.bg", input_filepath))))
#'     current_file = basename(gsub(input_extension, "__tempfile.bg", input_filepath))
#'   }
#'   
#'   # LiftOver bedGraph 
#'   system.time({system(paste(system.file("extdata", "liftOver", package = "genomeTools"), current_file, chain_file, basename(gsub(input_extension, "_liftover__tempfile.bg", input_filepath)), "/dev/null"))
#'   current_file = basename(gsub(input_extension, "_liftover__tempfile.bg", input_filepath))})
#'   
#'   # Remove duplicated regions post liftover and sort
#'   system(paste("sort -k1,1 -k2,2n -o", current_file, current_file))
#'   system(paste(system.file("extdata", "bedRemoveOverlap", package = "genomeTools"), current_file, gsub("_liftover__tempfile.bg", "liftover_duplicated_regions_removed__tempfile.bg", current_file)))
#'   current_file = gsub("_liftover__tempfile.bg", "liftover_duplicated_regions_removed__tempfile.bg", current_file)
#'   
#'   # If output_filetype is bg, can return current_file. If output_filetype is bw, need to convert it back to a bigwig
#'   if(output_filetype == "bg"){
#'     system(paste("mv", current_file, output_filename))
#'   } else if(output_filetype == "bw"){
#'     system(paste(system.file("extdata", "bedGraphToBigWig", package = "genomeTools"), current_file, chrom_sizes, output_filename))
#'   }
#'   
#'   # Remove temporary files
#'   system("rm *__tempfile*")
#'   
#' }
#'   
#'   
#'   
#'  
#' 
#' 
#' 
#' 
#' 
#' 
