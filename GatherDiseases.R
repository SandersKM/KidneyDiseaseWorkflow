start_time <- Sys.time()
library(rvest)
library(httr)
library(rentrez)
library(data.table)
library(XML)
library(xlsx)

##############
# User Input
##############

# Put Disease/Phenotype of interest below
disease.keyword = "Kidney Disease"

# Enter the path you would like the final CSV to be in:
disease_file_path = "/Users/ksanders/Documents/"


# Set to true what file formats you want:
write_to_excel <- TRUE
write_to_csv <- TRUE
write_to_tsv <- TRUE

disease_file_name = paste("Mendelian_Pipeline_",sub(" ", "_", disease.keyword), sep = "")

################
# MedGen Search
################

# API lookup
disease.search.count <- entrez_search(db="medgen", term=disease.keyword, retmax=0)$count
disease.search <- entrez_search(db="medgen", term=disease.keyword, retmax=disease.search.count)

# makes a data frame for the disease results
disease_file <- data.frame(id=disease.search$ids, stringsAsFactors = FALSE)
disease_file$summary <- entrez_summary(db="medgen", id =disease_file$id)
disease_file$conceptid <- extract_from_esummary(disease_file$summary, "conceptid")
disease_file$title <- extract_from_esummary(disease_file$summary, "title")
disease_file$definition <- sapply(extract_from_esummary(disease_file$summary, "definition"), '[', 1)
disease_file$definition <- unlist(lapply(disease_file$definition, toString))
disease_file <- disease_file[ , !(names(disease_file) %in% c("summary"))]
disease_genes <- entrez_link(dbfrom = "medgen", id = disease_file$id, db= 'gene')$links[1]$medgen_gene_diseases
disease_genes_symbols <- extract_from_esummary(entrez_summary(db="gene", id = disease_genes), c("uid", "name"))
disease_gene_by_id <- entrez_link(dbfrom = "medgen", id = disease_file$id, db= 'gene',by_id = TRUE)
disease_file$genes_symbol <- character(dim(disease_file)[1])
disease_file$genes <- character(dim(disease_file)[1])
# put together genes and gene symbol information for each disease 
for( i in 1:dim(disease_file)[1]){
  rownum <- match(xmlValue(disease_gene_by_id[i][[1]]$file[2]$IdList[1]$Id[1]$text), disease_file$id)
  disease_file$genes[rownum] <-
    paste0(disease_gene_by_id[i][[1]]$links$medgen_gene_diseases, collapse = "; ")
  disease_file$genes_symbol[rownum] <- paste0(disease_genes_symbols["name",disease_gene_by_id[i][[1]]$links$medgen_gene_diseases], collapse = "; ") 
}

# take out diseases without associated genes
disease_file <- disease_file[!disease_file$genes == "",] 
# Get rid of Summary Column
disease_file <- disease_file[ , !(names(disease_file) %in% c("summary"))] 

# Get inheritance type from disease summary
get_inheritance <- function(n){
  def <- strsplit(tolower(n), split = " ")
  if(length(def) > 0){
    if("autosomal" %in% def[[1]]){
      if("dominant" %in% def[[1]]){
        return("AD")
      }
      if("recessive" %in% def[[1]]){
        return("AR")
      }
      return("A")
    }
    if(("x"  %in% def[[1]]) || ("linked"  %in% def[[1]])|| ("x-linked"  %in% def[[1]])){
      return("X")
    }
  }
  return("")
}
disease_file$inheritance <- sapply(disease_file$definition, get_inheritance)
disease_file$inheritance <- unlist(disease_file$inheritance)

#################
# Write to file
#################

# Write results to csv
if(write_to_csv){
  write.csv(disease_file, file=paste(disease_file_path, disease_file_name, ".csv", sep=""), row.names = FALSE)
}
if(write_to_tsv){
  write.table(disease_file, file = paste(disease_file_path, disease_file_name, ".txt", sep=""), sep = "\t", row.names = FALSE)
}
if(write_to_excel){
  write.xlsx2(disease_file, file=paste(disease_file_path, disease_file_name, ".xlsx", sep=""), sheetName = "Disease", row.names = FALSE)
}


end_time <- Sys.time()