library(rvest)
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
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# get ensemblIDs
get_ensemblID <- function(n){
  id <- getBM("ensembl_gene_id", filters="hgnc_symbol", values=gene_file$name[n], mart=ensembl)
  return(id)
}
gene_file$EnsemblID <- sapply(1:dim(gene_file)[1], get_ensemblID)
# get start position of gene for hg19
get_start_position <- function(n){
  id <- getBM("start_position", filters="hgnc_symbol", values=gene_file$name[n], mart=ensembl)
  return(as.numeric(id[[1]]))
}
gene_file$start_position <- lapply(1:dim(gene_file)[1], get_start_position)
# get end position of gene for hg19
get_end_position <- function(n){
  id <- getBM("end_position", filters="hgnc_symbol", values=gene_file$name[n], mart=ensembl)
  return(as.numeric(id[[1]]))
}
gene_file$end_position <- sapply(1:dim(gene_file)[1], get_end_position)
# get percentage of GC content
get_percentage_gc_content <- function(n){
  id <- getBM("percentage_gene_gc_content", filters="hgnc_symbol", values=gene_file$name[n], mart=ensembl)
  return(id[[1]])
}
gene_file$percentage_gc_content <- sapply(1:dim(gene_file)[1], get_percentage_gc_content)

gene_file <- gene_file[c("name", "description", "phenotype", "summary", "geneID", "EnsemblID", "mim",
                 "map.location", "chromosome", "start_position", "end_position", "exon.count", "percentage_gc_content")]

# Web scrapping on Human Protein Atlas for protein/rna expression data
gene_file$hpa.url <- sapply(gene_file$EnsemblID, function(x){
  if(!is.logical(x)){
    # if there is > 1 EnsemblID, I am just using the first. 
    if(length(x) > 1){
      paste("https://www.proteinatlas.org/",x[[1]],"/tissue/kidney", sep = "")
    } 
    else{
      paste("https://www.proteinatlas.org/",x,"/tissue/kidney", sep = "")
    }
  }
}) 

gene_file$page <- sapply(gene_file$hpa.url, function(x){ 
  print(x)
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
gene_file$hpa.rna.gtex.rpkm <- sapply(gene_file$hpa.rna.expression, function(x){
  if(!is.null(x)){
    tryCatch({as.list(x[2] %>% html_nodes("td") %>% html_text())[[2]]},
             error = function(e){
               return(NULL)
             })}})
gene_file$hpa.rna.fantom5.tagspermillion <- sapply(gene_file$hpa.rna.expression, function(x){
  if(!is.null(x)){
    tryCatch({as.list(x[3] %>% html_nodes("td") %>% html_text())[[2]]},
             error = function(e){
               return(NULL)
             })}})

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
gene_file$hpa.protein.glomeruli <- sapply(gene_file$hpa.protein.expression, function(x){
  if(!is.null(x)){
    tryCatch({as.list(x[1] %>% html_nodes("td") %>% html_text)[[2]]},
             error = function(e){
               return(NULL)
             })}})
gene_file$hpa.protein.tubules <- sapply(gene_file$hpa.protein.expression, function(x){
  if(!is.null(x)){
    tryCatch({as.list(x[2] %>% html_nodes("td") %>% html_text)[[2]]},
             error = function(e){
               return(NULL)
             })}}) 
gene_file <- gene_file[ , !(names(gene_file) %in% c("page", "hpa.rna.expression",
                                                    "hpa.protein.expression"))]

# Web scrapping to get gene structure exons 
# let me know if you also want UTRs
gene_file$rcsb.url <- sapply(gene_file$name, function(x){
  paste("http://www.rcsb.org/pdb/gene/",x,"?v=hg19", sep = "")
}) 
gene_file$rcsb.genestructure.txt <- sapply(gene_file$rcsb.url, function(x){ 
  tryCatch({
    read_html(x) %>% 
      html_nodes("body div.container table.table.table-hover tr td") %>% 
      html_text()
  },
           error = function(e){
             return(NULL)})
  })
gene_file$exon <- sapply(1:dim(gene_file)[1], get_exons)

get_exons <- function(n){
  if(!(length(gene_file$exon.count[[n]]) == 0)){
    txt <- gene_file$rcsb.genestructure.txt[n]
    txt <- txt[[1]]
    exons <- vector(mode = "list", length = as.numeric(gene_file$exon.count[[n]]))
    i <- 1
    e <- 1
    while(i <= length(txt)){
      if("Exon" %in% txt[i]){
        exons[[e]] <- c(txt[i + 1], txt[i + 2])
        i <- i + 3
        e <- e + 1
      }
      else{
        i <- i + 1
      }
    }
    return(exons[!sapply(exons, is.null)])
  }
  else{return(NULL)}
  
}

gene_file <- gene_file[ , !(names(gene_file) %in% c("exon.count", "rcsb.genestructure.txt"))]

# for some reason, not all of the exon counts match the number of exon regions retrieved?
# this code is to make the two match...
gene_file$exon.count <- sapply(1:length(gene_file$exon)[1], function(n){
  if(!is.null(gene_file$exon[n])){
  length(gene_file$exon[[n]])}})



# UniProt/swissprot
# get_uniprotswissprot <- function(n){
#   id <- getBM("uniprotswissprot", filters="hgnc_symbol", values=gene_file$name[n], mart=ensembl)
#   return(id)
# }
# gene_file$uniprotswissprot<- sapply(1:dim(gene_file)[1], get_uniprotswissprot)

