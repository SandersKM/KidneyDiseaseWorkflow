
library(biomaRt)

# Go to https://www.ncbi.nlm.nih.gov/gene/advanced to search for all
# genes with a certain phenotype for humans. 
# ex: (Kidney Disease[Disease/Phenotype]) AND Homo Sapiens[Organism] 
# Click "Send to" in the top right corner, select "File" with 
# Format = "Tabular (text)". Click "Create File".

# Path to the downloaded file (include file name)
file_name <- "/Users/ksanders/Downloads/Kidney_Genes.txt"
gene_file <- read.delim(file_name, header = TRUE, sep = "\t")

# drops unneeded/repetitive columns
drops <- c("tax_id", "Org_name", "CurrentID", "X", "Aliases", "Status")
gene_file <- gene_file[ , !(names(gene_file) %in% drops)]

# set up biomart
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# get ensemblIDs
get_ensemblID <- function(n){
  id <- getBM("ensembl_gene_id", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  return(id)
}
gene_file$EnsemblID <- sapply(1:dim(gene_file)[1], get_ensemblID)
# get percentage of GC content
get_percentage_gc_content <- function(n){
  id <- getBM("percentage_gene_gc_content", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  return(id)
}
gene_file$percentage_gc_content <- sapply(1:dim(gene_file)[1], get_percentage_gc_content)
# get phenotype description
get_phenotype_description <- function(n){
  id <- getBM("phenotype_description", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  return(id)
}
gene_file$phenotype_description <- sapply(1:dim(gene_file)[1], get_phenotype_description)

# MIM disease description

get_MIM_disease_description <- function(n){
  id <- getBM("mim_morbid_description", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  return(id)
}
gene_file$MIM_disease<- sapply(1:dim(gene_file)[1], get_MIM_disease_description)

# UniProt/swissprot

get_uniprotswissprot <- function(n){
  id <- getBM("uniprotswissprot", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  return(id)
}
gene_file$uniprotswissprot<- sapply(1:dim(gene_file)[1], get_uniprotswissprot)


# WikiGene Description

get_wikigene_description <- function(n){
  id <- getBM("wikigene_description", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  return(id)
}
gene_file$wikigene_description <- sapply(1:dim(gene_file)[1], get_wikigene_description)





# # get go ID
# get_go_id <- function(n){
#   id <- getBM("go_id", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
#   return(id)
# }
# gene_file$go_id <- sapply(1:dim(gene_file)[1], get_go_id)
# # get go name
# get_name_1006 <- function(n){
#   id <- getBM("name_1006", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
#   return(id)
# }
# gene_file$name_1006 <- sapply(1:dim(gene_file)[1], get_name_1006)
# # get go definition
# get_definition_1006 <- function(n){
#   id <- getBM("definition_1006", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
#   return(id)
# }
# gene_file$definition_1006 <- sapply(1:dim(gene_file)[1], get_definition_1006)

# get GOSlim GOA Accession(s)
# get_goslim_goa_accession <- function(n){
#   id <- getBM("goslim_goa_accession", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
#   return(id)
# }
# gene_file$goslim_goa_accession <- sapply(1:dim(gene_file)[1], get_goslim_goa_accession)
# # get GOSlim GOA Descriptions
# get_goslim_goa_description <- function(n){
#   id <- getBM("goslim_goa_description", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
#   return(id)
# }
# gene_file$goslim_goa_description <- sapply(1:dim(gene_file)[1], get_goslim_goa_description)
