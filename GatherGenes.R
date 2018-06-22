
library(biomaRt)
library(gwascat)
data(ebicat37)
library(rentrez)

gene_file <- data.frame(geneID = integer(dim(disease_file)[1]), phenotype = character(dim(disease_file)[1]), 
                        stringsAsFactors = FALSE)
index <- 1
for(i in 1:dim(disease_file)[1]){
  spl <- strsplit(disease_file$genes[i], split = "; ")
  firstid <- as.integer(spl[[1]][1])
  if(firstid %in% gene_file$geneID){
    idrow <- which(disease_file$genes == firstid, arr.ind = TRUE)
    gene_file[idrow, 2] = paste(gene_file[idrow, 2],disease_file$title[i], sep = "; ")
  }
  else{
    gene_file[index, 1] = firstid
    gene_file[index, 2] = as.character(disease_file$title[i])
    index <- index + 1
  } 
  if(length(spl[[1]]) > 1){
    for(j in 2:length(spl[[1]])){
      if(index < i){
        gene_file[index, 1] = as.integer(spl[[1]][j])
        gene_file[index, 2] =as.character(disease_file$title[i])
        index <- index + 1
      }
      else{
        tempdf <- data.frame(geneID = as.integer(spl[[1]][j]), phenotype = disease_file$title[i])
        gene_file <- rbind(gene_file, tempdf)
      }
    }
  }
}
gene_file <- gene_file[!gene_file$geneID == 0,]
rm(tempdf)










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
# get start position of gene for hg19
get_start_position <- function(n){
  id <- getBM("start_position", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  return(as.numeric(id))
}
gene_file$start_position_on_the_genomic_accession <- sapply(1:dim(gene_file)[1], get_start_position)
# get end position of gene for hg19
get_end_position <- function(n){
  id <- getBM("end_position", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  return(as.numeric(id))
}
gene_file$end_position_on_the_genomic_accession <- sapply(1:dim(gene_file)[1], get_end_position)
# get percentage of GC content
get_percentage_gc_content <- function(n){
  id <- getBM("percentage_gene_gc_content", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  return(id)
}
gene_file$percentage_gc_content <- sapply(1:dim(gene_file)[1], get_percentage_gc_content)

# get phenotype description
get_phenotype_description <- function(n){
  id <- getBM("phenotype_description", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  if(length(id) == 0){return("NA")}
  return(id)
}
gene_file$phenotype_description <- sapply(1:dim(gene_file)[1], get_phenotype_description)
gene_file$phenotype_description <- sapply(gene_file$phenotype_description, function(x) {paste(x, collapse = "; ")})



# UniProt/swissprot

get_uniprotswissprot <- function(n){
  id <- getBM("uniprotswissprot", filters="hgnc_symbol", values=gene_file$Symbol[n], mart=ensembl)
  return(id)
}
gene_file$uniprotswissprot<- sapply(1:dim(gene_file)[1], get_uniprotswissprot)





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

gwasgene <- subsetByTraits(ebicat37, tr="Chronic kidney disease"  )
