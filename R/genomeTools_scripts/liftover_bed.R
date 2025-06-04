#' #' Liftover a bed file from one genome build to another
#' #'
#' #' Liftover a bed file from one genome build to another.
#' #' bedRemoveOverlap is used to remove overlapping regions resulting from the output files after liftover. 
#' #' Output files have the basename of the input files with "liftover" and the name of the target geneome build (e.g. "hg38") depending on value of to_genome_build 
#' #'
#' #' @param input_filepaths Paths to BED files to be lifted-over. BED files may be zipped.
#' #' @param to_genome_build Either "hg18", hg19" or "hg38" depending on the intended output genome.
#' #' @param output_directory Directory to save output files in. Default is current working directory. If output_directory doesn't exist, it is created.
#' #' @param compress_output_beds A logical value indicating whether output bedGraph files should be compressed. Default value is True
#' #' @return None
#' #' @export
#' liftover_bed = function(input_filepaths, to_genome_build, output_directory = ".", compress_output_beds = T){
#'   
#'   # Check that allowed values are entered for input_filetype, output_filetype, to_genome_build and ncores
#'   match.arg(arg = to_genome_build, choices = c("hg18", "hg19", "hg38"), several.ok = F)
#'   
#'   # Create output_directory if it doesn't exist)
#'   dir.create(output_directory, recursive = T, showWarnings = F)
#'   
#'   `%dopar%` = foreach::`%dopar%`
#'   
#'   # Create output filenames 
#'   output_filenames = paste(output_directory, paste0(tools::file_path_sans_ext(gsub(".gz", "", input_filepaths)), "_liftover_", to_genome_build, ".", tools::file_ext(input_filepaths)), sep = "/")
#'   
#'   # Load correct chain file for relevent genome build
#'   if (to_genome_build == "hg18"){
#'     chain_file = system.file("extdata", "hg18ToHg19.over.chain", package = "genomeTools")
#'     chrom_sizes = system.file("extdata", "hg19_chrom_size.tsv", package = "genomeTools")
#'   } else if (to_genome_build == "hg19"){
#'     chain_file = system.file("extdata", "hg38ToHg19.over.chain", package = "genomeTools")
#'     chrom_sizes = system.file("extdata", "hg19_chrom_size.tsv", package = "genomeTools")
#'   } else if (to_genome_build == "hg38"){
#'     chain_file = system.file("extdata", "hg19ToHg38.over.chain", package = "genomeTools")
#'     chrom_sizes = system.file("extdata", "hg38_chrom_size.tsv", package = "genomeTools")
#'   }
#'   
#'   # Liftover files
#'   foreach::foreach(file_number = seq_along(input_filepaths)) %dopar% {
#'   
#'     # Set current working file to the input file
#'     current_file = input_filepaths[file_number]
#'     
#'     # LiftOver bed
#'     liftover_temp = tempfile()
#'     system(paste(system.file("extdata", "liftOver", package = "genomeTools"), current_file, chain_file, liftover_temp, "/dev/null"))
#'     current_file = liftover_temp
#'     
#'     # Remove duplicated regions post liftover and sort
#'     system("echo Removing duplicated regions")
#'     overlap_removed_temp = tempfile()
#'     system(paste("sort -k1,1 -k2,2n -o", current_file, current_file))
#'     system(paste(system.file("extdata", "bedRemoveOverlap", package = "genomeTools"), current_file, overlap_removed_temp))
#'     current_file = overlap_removed_temp
#'     unlink(liftover_temp)liftover_bed = function(input_filepaths, to_genome_build, output_directory = ".", compress_output_beds = T){
#'   
#'   # Check that allowed values are entered for input_filetype, output_filetype, to_genome_build and ncores
#'   match.arg(arg = to_genome_build, choices = c("hg19", "hg38"), several.ok = F)
#'   
#'   # Create output_directory if it doesn't exist)
#'   dir.create(output_directory, recursive = T, showWarnings = F)
#'   
#'   `%dopar%` = foreach::`%dopar%`
#'   
#'   # Create output filenames 
#'   output_filenames = paste(output_directory, paste0(tools::file_path_sans_ext(gsub(".gz", "", input_filepaths)), "_liftover_", to_genome_build, ".", tools::file_ext(input_filepaths)), sep = "/")
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
#'   # Liftover files
#'   foreach::foreach(file_number = seq_along(input_filepaths)) %dopar% {
#'   
#'     # Set current working file to the input file
#'     current_file = input_filepaths[file_number]
#'     
#'     # LiftOver bed
#'     liftover_temp = tempfile()
#'     system(paste(system.file("extdata", "liftOver", package = "genomeTools"), current_file, chain_file, liftover_temp, "/dev/null"))
#'     current_file = liftover_temp
#'     
#'     # Remove duplicated regions post liftover and sort
#'     system("echo Removing duplicated regions")
#'     overlap_removed_temp = tempfile()
#'     system(paste("sort -k1,1 -k2,2n -o", current_file, current_file))
#'     system(paste(system.file("extdata", "bedRemoveOverlap", package = "genomeTools"), current_file, overlap_removed_temp))
#'     current_file = overlap_removed_temp
#'     unlink(liftover_temp)
#'     
#'     # Move file to specified destination
#'     system(paste("mv", current_file, output_filenames[file_number]))
#'     
#'     # Gzip output BED if specified
#'     if(compress_output_beds){
#'       system(paste("gzip", output_filenames))
#'     }
#'     
#'     # Remove temporary files
#'     unlink(overlap_removed_temp)
#'   }
#'   
#' }
#'     
#'     # Move file to specified destination
#'     system(paste("mv", current_file, output_filenames[file_number]))
#'     
#'     # Gzip output BED if specified
#'     if(compress_output_beds){
#'       system(paste("gzip", output_filenames))
#'     }
#'     
#'     # Remove temporary files
#'     unlink(overlap_removed_temp)
#'   }
#'   
#' }