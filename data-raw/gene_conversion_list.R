### Make a list which can be used to convert between gene IDs using bioMart

# Load required packages
library(biomaRt)

# Get gene mart for hg38 version 103 
gene_mart_hg38 = useEnsembl("ENSEMBL", dataset="hsapiens_gene_ensembl", version = 103)

# Download a data.frame with ensembl gene IDs, HGNC IDs and entrezgene IDs. Downloaded 2/4/21
system.time({gene_match_df = getBM(mart = gene_mart_hg38, 
    attributes = c("ensembl_gene_id", "hgnc_id", "hgnc_symbol", "entrezgene_id"))})

# Make entrezgene_id a character vector
gene_match_df$entrezgene_id = as.character(gene_match_df$entrezgene_id)

# Remove HGNC: from start of HGNC IDs
gene_match_df$hgnc_id = gsub("HGNC:", "", gene_match_df$hgnc_id)

# Set HGNC IDs and symbols with a value of "" = to NA
gene_match_df$hgnc_id[!nzchar(gene_match_df$hgnc_id)] = NA
gene_match_df$hgnc_symbol[!nzchar(gene_match_df$hgnc_symbol)] = NA

# Remove rows where three of the entries have NA values
gene_match_df = gene_match_df[which(apply(gene_match_df, 1, function(x) sum(is.na(x))) < 3), ]

# Create a list of names vectors matching different gene ID combinations. If there are multiple matches between identifiers, select the first match in lexicographical order
gene_conversion_list  = list(
  ensembl_id_hgnc_id = sapply(split(gene_match_df$hgnc_id, gene_match_df$ensembl_gene_id), function(x) sort(x)[1]),
  ensembl_id_hgnc_symbol = sapply(split(gene_match_df$hgnc_symbol, gene_match_df$ensembl_gene_id), function(x) sort(x)[1]),
  ensembl_id_entrez_id = sapply(split(gene_match_df$entrezgene_id, gene_match_df$ensembl_gene_id), function(x) sort(x)[1]),
  hgnc_id_ensembl_id = sapply(split(gene_match_df$ensembl_gene_id, gene_match_df$hgnc_id), function(x) sort(x)[1]),
  hgnc_id_hgnc_symbol = sapply(split(gene_match_df$hgnc_symbol, gene_match_df$hgnc_id), function(x) sort(x)[1]),
  hgnc_id_entrez_id = sapply(split(gene_match_df$entrezgene_id, gene_match_df$hgnc_id), function(x) sort(x)[1]),
  hgnc_symbol_ensembl_id = sapply(split(gene_match_df$ensembl_gene_id, gene_match_df$hgnc_symbol), function(x) sort(x)[1]),
  hgnc_symbol_hgnc_id = sapply(split(gene_match_df$hgnc_id, gene_match_df$hgnc_symbol), function(x) sort(x)[1]),
  hgnc_symbol_entrez_id = sapply(split(gene_match_df$entrezgene_id, gene_match_df$hgnc_symbol), function(x) sort(x)[1]),
  entrez_id_ensembl_id = sapply(split(gene_match_df$ensembl_gene_id, gene_match_df$entrezgene_id), function(x) sort(x)[1]),
  entrez_id_hgnc_id = sapply(split(gene_match_df$hgnc_id, gene_match_df$entrezgene_id), function(x) sort(x)[1]),
  entrez_id_hgnc_symbol = sapply(split(gene_match_df$hgnc_symbol, gene_match_df$entrezgene_id), function(x) sort(x)[1])
)
save(gene_conversion_list, file = "~/my_packages/genes/data/gene_conversion_list.Rdata")