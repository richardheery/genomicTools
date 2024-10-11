#' Convert gene identifiers
#'
#' Convert between Ensembl IDs, HGNC IDs, HGNC symbols and Entrez IDs. 
#' Matches IDs using Ensembl version 112. 
#'
#' @param gene_ids A character vector with gene IDs to convert
#' @param input_id_type Type of input gene ID. Options are "ensembl_id", "hgnc_id", "hgnc_symbol" or "entrez_id". 
#' @param output_id_type Type of output gene ID. Options are "ensembl_id", "hgnc_id", "hgnc_symbol" or "entrez_id".
#' @return None
#' @export
gene_id_convert = function(gene_ids, input_id_type, output_id_type){
  
  # Check that allowed choices are selected for input ID and output ID and that they are different
  match.arg(arg = input_id_type, choices = c("ensembl_id", "hgnc_id", "hgnc_symbol", "entrez_id"), several.ok = F)
  match.arg(arg = output_id_type, choices = c("ensembl_id", "hgnc_id", "hgnc_symbol", "entrez_id"), several.ok = F)
  if(input_id_type == output_id_type){stop("Input ID and output ID type must be different")}
  
  # Ensure gene IDs are a character vector
  gene_ids = as.character(gene_ids)
  
  # If input ID type is Ensembl IDs, remove the version number if present
  if (input_id_type == "ensembl_id"){gene_ids = gsub("\\..*", "", gene_ids)}
  
  # Load gene_conversion_list
  gene_conversion_list = genomicTools:::gene_conversion_list
  
  # Convert genes and return
  if (input_id_type == "ensembl_id" && output_id_type == "hgnc_id") {return(gene_conversion_list$ensembl_id_hgnc_id[gene_ids])}
  else if (input_id_type == "ensembl_id" && output_id_type == "hgnc_symbol") {return(gene_conversion_list$ensembl_id_hgnc_symbol[gene_ids])}
  else if (input_id_type == "ensembl_id" && output_id_type == "entrez_id") {return(gene_conversion_list$ensembl_id_entrez_id[gene_ids])}
  else if (input_id_type == "hgnc_id" && output_id_type == "ensembl_id") {return(gene_conversion_list$hgnc_id_ensembl_id[gene_ids])}
  else if (input_id_type == "hgnc_id" && output_id_type == "hgnc_symbol") {return(gene_conversion_list$hgnc_id_hgnc_symbol[gene_ids])}
  else if (input_id_type == "hgnc_id" && output_id_type == "entrez_id") {return(gene_conversion_list$hgnc_id_entrez_id[gene_ids])}
  else if (input_id_type == "hgnc_symbol" && output_id_type == "ensembl_id") {return(gene_conversion_list$hgnc_symbol_ensembl_id[gene_ids])}
  else if (input_id_type == "hgnc_symbol" && output_id_type == "hgnc_id") {return(gene_conversion_list$hgnc_symbol_hgnc_id[gene_ids])}
  else if (input_id_type == "hgnc_symbol" && output_id_type == "entrez_id") {return(gene_conversion_list$hgnc_symbol_entrez_id[gene_ids])}
  else if (input_id_type == "entrez_id" && output_id_type == "ensembl_id") {return(gene_conversion_list$entrez_id_ensembl_id[gene_ids])}
  else if (input_id_type == "entrez_id" && output_id_type == "hgnc_id") {return(gene_conversion_list$entrez_id_hgnc_id[gene_ids])}
  else if (input_id_type == "entrez_id" && output_id_type == "hgnc_symbol") {return(gene_conversion_list$entrez_id_hgnc_symbol[gene_ids])}
  
}