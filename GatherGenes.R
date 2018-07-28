start_time <- Sys.time()
library(rvest)
library(biomaRt)
library(rentrez)
suppressPackageStartupMessages(library(GenomicFeatures))
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Enter the path you would like the final CSV to be in:
gene_file_path = "/Users/ksanders/Documents/"


gene_file <- data.frame(geneID = unique(unlist(strsplit(disease_file$genes, split = "; "))), 
                        phenotype = character(dim(disease_file)[1]), 
                        stringsAsFactors = FALSE)

split_gene <- strsplit(disease_file$genes, split = "; ")
for(i in 1:length(split_gene)){
  for( j in 1:length(split_gene[[i]])){
    rownum <- which(gene_file$geneID == split_gene[i][[1]][j])
    if(gene_file$phenotype[rownum] !=""){
      gene_file$phenotype[rownum] <- 
        paste(gene_file$phenotype[rownum],disease_file$title[i],sep = "; ")
    }
    else{
      gene_file$phenotype[rownum] <- disease_file$title[i]
    }
  }
}
gene_file <- gene_file[!gene_file$geneID == 0,]


gene_file$summary <- entrez_summary(db="gene", id =gene_file$geneID)
gene_file$name <- extract_from_esummary(gene_file$summary, "name")
gene_file$description <- extract_from_esummary(gene_file$summary, "description")
gene_file$chromosome <- extract_from_esummary(gene_file$summary, "chromosome")
gene_file$map.location <- extract_from_esummary(gene_file$summary, "maplocation")
gene_file$mim <- extract_from_esummary(gene_file$summary, "mim")
gene_file$exon.count <- sapply(1:dim(gene_file)[1], function(x){as.integer(extract_from_esummary(
  gene_file$summary,"genomicinfo")[[x]]$exoncount)})
gene_file$exon.count <- sapply(gene_file$exon.count, as.numeric)
gene_file$summary <- extract_from_esummary(gene_file$summary, "summary")

# set up biomart
if(!exists("ensembl")){
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
}
# get information from biomart
allBM <- getBM(c("hgnc_symbol","ensembl_gene_id","start_position", "end_position", "strand"), filters="hgnc_symbol", 
                  values=gene_file$name, mart=ensembl)
# merge rows where 1 hgnc symbol has >1 ensembl ID
allBM <- aggregate(allBM[,-1], list(allBM[,1]), function(x) paste0(unique(x), collapse = "; "))
names(allBM)[names(allBM) == 'Group.1'] <- "name"
# merge the biomart results with the gene file
gene_file = merge(x = gene_file, y = allBM, by="name",all.x=T, all.y = F)


# Web scrapping on Human Protein Atlas for protein/rna expression data in kidney
# If you are not looking at kidney disease, this can be easily modified. 
gene_file$hpa.url <- sapply(gene_file$ensembl_gene_id, function(x){
  if(!is.logical(x)){
    # if there is > 1 ensembl_gene_id, I am just using the first. 
    if(length(x) > 1){
      paste("https://www.proteinatlas.org/",x[[1]],"/tissue/kidney", sep = "")
    } 
    else{
      paste("https://www.proteinatlas.org/",x,"/tissue/kidney", sep = "")
    }
  }
}) 

gene_file$page <- sapply(gene_file$hpa.url, function(x){ 
  if(!is.null(x)){
    tryCatch({read_html(x)},
             error = function(e){
               return(NULL)
             })
  }
})

gene_file$hpa.rna.expression <- sapply(gene_file$page, function(x){
  if(!is.null(x)){
    x %>% html_nodes("body table.main_table tr div.menu_margin 
                     table.border.dark.round table.noborder.nowrap tr"
                     )}})

gene_file$hpa.rna.hpa.tpm <- sapply(gene_file$hpa.rna.expression, function(x){
  if(!is.null(x)){
    as.list(x[1] %>% html_nodes("td") %>% html_text())[[2]]}})
gene_file$hpa.rna.hpa.tpm <- lapply(gene_file$hpa.rna.hpa.tpm, toString)
gene_file$hpa.rna.hpa.tpm <- unlist(gene_file$hpa.rna.hpa.tpm)

gene_file$hpa.rna.gtex.rpkm <- sapply(gene_file$hpa.rna.expression, function(x){
  if(!is.null(x)){
    tryCatch({as.list(x[2] %>% html_nodes("td") %>% html_text())[[2]]},
             error = function(e){
               return(NULL)
             })}})
gene_file$hpa.rna.gtex.rpkm <- lapply(gene_file$hpa.rna.gtex.rpkm, toString)
gene_file$hpa.rna.gtex.rpkm <- unlist(gene_file$hpa.rna.gtex.rpkm)

gene_file$hpa.rna.fantom5.tagspermillion <- sapply(gene_file$hpa.rna.expression, function(x){
  if(!is.null(x)){
    tryCatch({as.list(x[3] %>% html_nodes("td") %>% html_text())[[2]]},
             error = function(e){
               return(NULL)
             })}})
gene_file$hpa.rna.fantom5.tagspermillion <- lapply(gene_file$hpa.rna.fantom5.tagspermillion, toString)
gene_file$hpa.rna.fantom5.tagspermillion <- unlist(gene_file$hpa.rna.fantom5.tagspermillion)

gene_file$hpa.protein.expression <- sapply(gene_file$page, function(x){
  if(!is.null(x)){
    x %>% html_nodes("body table.main_table tr div.menu_margin 
                     table table.dark th.nopadd table.border.dark table.noborder.nowrap tr")
    }})

gene_file$hpa.protein.glomeruli <- sapply(gene_file$hpa.protein.expression, function(x){
  if(!is.null(x)){
    tryCatch({as.list(x[1] %>% html_nodes("td") %>% html_text)[[2]]},
             error = function(e){
               return(NULL)
             })}})
gene_file$hpa.protein.glomeruli <- lapply(gene_file$hpa.protein.glomeruli, toString)
gene_file$hpa.protein.glomeruli <- unlist(gene_file$hpa.protein.glomeruli)

gene_file$hpa.protein.tubules <- sapply(gene_file$hpa.protein.expression, function(x){
  if(!is.null(x)){
    tryCatch({as.list(x[2] %>% html_nodes("td") %>% html_text)[[2]]},
             error = function(e){
               return(NULL)
             })}}) 
gene_file$hpa.protein.tubules <- lapply(gene_file$hpa.protein.tubules, toString)
gene_file$hpa.protein.tubules <- unlist(gene_file$hpa.protein.tubules)
gene_file <- gene_file[ , !(names(gene_file) %in% c("page", "hpa.rna.expression",
                                                    "hpa.protein.expression"))]

gene_file$cannonical_transcript <- character(dim(gene_file)[1])
gene_file$exon <- character(dim(gene_file)[1])
get_exons <- function(n){
  print(gene_file$name[n])
  try({
    transcripts <- as.data.frame(fromJSON(paste("http://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/", gene_file$name[n],
                                                "?content-type=application/json;expand=1", sep = "")))
    canonical.rownum <- which(transcripts$Transcript.is_canonical == 1)
    gene_file$cannonical_transcript[n] <<- transcripts$Transcript.id[canonical.rownum]
    gene_file$exon[n] <<- paste0(transcripts$Transcript.Exon[canonical.rownum][[1]]$start, " : ",
                             transcripts$Transcript.Exon[canonical.rownum][[1]]$end, collapse = "; ")
    })
}

sapply(1:dim(gene_file)[1], get_exons)

get_gnomad_website_gene <- function(n){
  if(class(gene_file$ensembl_gene_id[[n]]) != "logical"){
    base_url <- "http://gnomad.broadinstitute.org/gene/"
    variantID <- gene_file$ensembl_gene_id[n]
    return(paste(base_url,variantID,sep = ""))
  }
  return(NULL)
}
gene_file$gnomAD.website <- sapply(1:dim(gene_file)[1], get_gnomad_website_gene)
gene_file$gnomAD.website <- lapply(gene_file$gnomAD.website, toString)
gene_file$gnomAD.website <- unlist(gene_file$gnomAD.website)

###############################################
# PubMed IDs containing Gene Symbol AND disease keywords
###############################################

# gene_file$pubids <- sapply(gene_file$name, function(x){
#   gsub(",",";",toString((entrez_search(db = "pubmed",
#                                        term = paste("(",x ," AND ", disease.keyword,")", 
#                                                    sep = ""), retmax = 9999)$ids)))})

###############################################
# Clinvar IDs for Pathogenic variants
###############################################

# gene_file$clinvar.pathogenic.variants <- sapply(gene_file$geneID, function(x){
#   gsub(",",";",toString(entrez_link(dbfrom = "gene", id = x, 
#                                     db="all")$links$gene_clinvar_specific))
# })

#################
# Write to file
#################

if(!exists(disease_file_path)){
  if(!exists(disease.keyword)){
    disease.keyword <- "?"
  }
  disease_file_path <-  paste("Mendelian_Pipeline_",sub(" ", "_", disease.keyword), sep = "")
}

# Write results to csv
if(write_to_csv){
  write.csv(gene_file, file=paste(disease_file_path, disease_file_name, "_genes", ".csv", sep=""), row.names = FALSE)
}
if(write_to_tsv){
  write.table(gene_file, file = paste(disease_file_path, disease_file_name, "_genes",".txt", sep=""), sep = "\t", row.names = FALSE)
}
if(write_to_excel){
  write.xlsx2(gene_file, file=paste(disease_file_path, disease_file_name, ".xlsx", sep=""), sheetName = "Genes", row.names = FALSE)
}
k<- disease_file_path

end_time <- Sys.time()