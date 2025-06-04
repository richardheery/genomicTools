#' Convert gene identifiers
#'
#' Convert between Ensembl IDs, HGNC IDs, HGNC symbols and Entrez IDs
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
  gene_conversion_list = genomeTools:::gene_conversion_list
  
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

#' Update gene symbols
#'
#' Update symbols using HGNChelper
#'
#' @param gene_symbols A vector of gene names to update
#' @param ambiguous_as_na A logical value indicating whether to replace ambiguous results (i.e. results mapping to more than one gene symbol) with NA. Default is FALSE.
#' @param remove_na A logical value indicating whether to remove NA values  and sort the returned results, Default if FALSE.
#' @return None
#' @export
update_gene_symbols = function(gene_symbols, ambiguous_as_na = F, remove_na = F){
  results = HGNChelper::checkGeneSymbols(gene_symbols, map = genomeTools:::hgnchelper_map_01_09_2020)$Suggested.Symbol
  if(ambiguous_as_na){results[grepl("///", results)] = NA}
  if(remove_na){results = sort(results)}
  return(results)
}

#' Selects the correct gene using chromosome name when multiple genes are returned by HGNChelper::checkGeneSymbols()
#'
#' Selects the correct gene for a given chromosome when multiple genes are returned by HGNChelper::checkGeneSymbols(), 
#' if one and only one of the returned genes is on the indicated chromosome. Returns NA if more than one or none of the genes are on the indicated chromosome. 
#'
#' @param gene_symbols A vector of gene names 
#' @param chromosomes A vector of the names of chromosomes associated with each gene in \code{gene_symbols}
#' @return None
#' @export
gene_select_by_chrom = function(genes, chromosomes){
  
  # Get list of chromosomes and the genes located on them from HGNC
  chromosome_gene_list = readRDS("~/genes/hgnc/chromosome_hgnc_genes_list.rds")
  
  # Find which of the returned genes are on the correct chromosome
  correct_genes = sapply(seq_along(genes), function(x) c(intersect(trimws(unlist(strsplit(genes[x], "///"))), 
    chromosome_gene_list[[chromosomes[x]]]), NULL))
  
  # Remove any results where more than one of the returned genes was on the correct chromosome or else no genes were on the correct chromosome
  correct_genes[lengths(correct_genes) != 1] = NA
  return(unlist(correct_genes))
}
