start_time <- Sys.time()
library(rvest)
library(httr)
library(rentrez)
library(data.table)

##########################################
# Put Disease/Phenotype of interest below
##########################################

disease.keyword = "Kidney Disease"

# Enter the path you would like the final CSV to be in:
disease_file_path = "/Users/ksanders/Documents/"
disease_file_name = paste(sub(" ", "_",disease.keyword),  ".csv", sep = "")


# API lookup
disease.search.count <- entrez_search(db="medgen", term=disease.keyword, retmax=0)$count
disease.search <- entrez_search(db="medgen", term=disease.keyword, retmax=disease.search.count)

# makes a data frame for the disease results
disease_file <- data.frame(id=disease.search$ids)
disease_file$summary <- entrez_summary(db="medgen", id =disease_file$id)
disease_file$conceptid <- extract_from_esummary(disease_file$summary, "conceptid")
disease_file$title <- extract_from_esummary(disease_file$summary, "title")
disease_file$definition <- sapply(extract_from_esummary(disease_file$summary, "definition"), '[', 1)
disease_file$definition <- unlist(lapply(disease_file$definition, toString))
disease_file <- disease_file[ , !(names(disease_file) %in% c("summary"))]
disease_genes <- entrez_link(dbfrom = "medgen", id = disease_file$id, db= 'gene')$links[1]$medgen_gene_diseases
disease_genes_symbols <- extract_from_esummary(entrez_summary(db="gene", id = disease_genes), c("uid", "name"))
demo_gene <- entrez_link(dbfrom = "medgen", id = disease_file$id, db= 'gene',by_id = TRUE)
disease_file$demo_gene <- character(dim(disease_file)[1])
get_genes <- function(n){
  return(entrez_link(dbfrom="medgen", id = disease_file$id[n], db='gene')$links)
}
disease_file$genes <- sapply(1:dim(disease_file)[1], get_genes)
disease_file$genes_symbol <- character(dim(disease_file)[1])
make_string_list <- function(n){
  if(length(disease_file$genes[[n]]) == 0){
    return("NULL")
  }
  disease_file$genes_symbol[n] <<- paste0(unlist(disease_genes_symbols["name",disease_file$genes[[n]]$medgen_gene_diseases]), collapse = "; ") 
  return(paste(disease_file$genes[[n]]$medgen_gene_diseases, collapse="; "))
}
# take out diseases without associated genes
disease_file$genes <- sapply(1:dim(disease_file)[1], make_string_list)
disease_file <- disease_file[!disease_file$genes == "NULL",]
disease_file <- disease_file[ , !(names(disease_file) %in% c("summary"))]

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

write.csv(disease_file, file=paste(disease_file_path, disease_file_name, sep=""), row.names = FALSE)

# you should go back and add in the names of the related genes. 
# Also make one excel file with multiple sheets
end_time <- Sys.time()