library(rvest)
library(httr)
library(rentrez)

########################################
# If you have an API key for NCBI
# talk about the api key setup and
# actually do it yourself ;)
#######################################



##########################################
# Put Disease/Phenotype of interest below
##########################################

disease.keyword = "Kidney Disease"

# does the API lookup
disease.search.count <- entrez_search(db="medgen", term=disease.keyword, retmax=0)$count
disease.search <- entrez_search(db="medgen", term=disease.keyword, retmax=disease.search.count)

# makes a data frame for the disease results
disease_file <- data.frame(id=disease.search$ids)
disease_file$summary <- entrez_summary(db="medgen", id =disease_file$id)
disease_file$conceptid <- extract_from_esummary(disease_file$summary, "conceptid")
disease_file$title <- extract_from_esummary(disease_file$summary, "title")
disease_file$definition <- sapply(extract_from_esummary(disease_file$summary, "definition"), '[', 1)
disease_file <- disease_file[ , !(names(disease_file) %in% c("summary"))]
get_genes <- function(n){
  links <- entrez_link(dbfrom="medgen", id = disease_file$id[n], db='gene')$links
  return(links)
}
disease_file$genes <- sapply(1:dim(disease_file)[1], get_genes)
make_string_list <- function(n){
  if(length(disease_file$genes[[n]]) == 0){
    return("NULL")
  }
  return(paste(disease_file$genes[[n]]$medgen_gene_diseases, collapse="; "))
}
# take out diseases without associated genes
disease_file$genes <- sapply(1:dim(disease_file)[1], make_string_list)
disease_file <- disease_file[!disease_file$genes == "NULL",]
disease_file <- disease_file[ , !(names(disease_file) %in% c("summary"))]

# Webscraping for inheritance and disease???
# gtr.gene.page <- read_html(paste("https://www.ncbi.nlm.nih.gov/gtr/genes/", gene_file$GeneID[3], "/", sep = ""))
# 
# gtr.disease.url <- html_nodes(gtr.gene.page, "td a")[1] %>% html_attr("href")
# gtr.disease.page <- read_html(paste("https://www.ncbi.nlm.nih.gov", gtr.disease.url, sep = ""))
# gtr.disease.name <- gtr.disease.page %>% html_node("h1") %>% html_text()
# gtr.disease.inheritance <- gtr.disease.page %>% html_node("div.grid div.col div.wrap div.page div.container div#maincontent 
#                                                           div.col1 div#gtr_page_cont div#gtr_maincontent div.rprt div.page_header dl dl dd a") %>% html_text()
# gtr.disease.inheritance <- html_nodes(gtr.disease.page, "dd a")[1] %>% html_text()
# gtr.disease.summary.long <- gtr.disease.page %>% html_nodes("div.grid div.col div.wrap div.page div.container div#maincontent 
#                                                             div.col1 div#gtr_page_cont div#gtr_maincontent div.rprt div.rprt-section div.rprt-section-body") %>% html_text()
# gtr.disease.summary <- gtr.disease.summary.long[1]