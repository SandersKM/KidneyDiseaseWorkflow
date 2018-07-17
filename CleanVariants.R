# Uncomment and run the following lines.
# It will take a few minutes, so be patient
# source("https://bioconductor.org/biocLite.R")
# biocLite("fitCons.UCSC.hg19")
# biocLite("phastCons100way.UCSC.hg19")

suppressPackageStartupMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(GenomicScores))
suppressPackageStartupMessages(library(phastCons100way.UCSC.hg19))
suppressPackageStartupMessages(library(fitCons.UCSC.hg19))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))

library(RISmed)
library(gwascat)
library(jsonlite)
library(rentrez)
library(httr)
library(xml2)
library(GenomeGraphs)
library(Gviz)

# download the variant CSV from the gnomAD browser for the 
# gene of interest. Write the path file below
all_variants_path <- "/Users/ksanders/Downloads/gnomad_ENSG00000185499_2018_07_10_08_45_10.csv"
variants <- read.csv(all_variants_path, header = TRUE, sep=",")
# Change the symbol to correspond to the gene of interest
gene.of.interest.symbol <- "MUC1"
# Enter the path you would like the final CSV to be in:
variants_file_path = "/Users/ksanders/Documents/"
variants_file_name = paste("variants_", gene.of.interest.symbol, ".csv", sep = "")





# function to get rid of variants flagged with "LC LoF" or "SEGDUP" or with low and modifier annotations
low_annotations <- c("splice region","synonymous","3' UTR","5' UTR","downstream gene","intron","upstream gene",
                     "non coding transcript exon")

clean <- function(n){
  if(variants$Annotation[n] %in% low_annotations){
    return(FALSE)
  }
  if(variants$Flags[n] != ""){
    return(FALSE)
  }
  return(TRUE)
}
variants$clean <- sapply(1:dim(variants)[1], clean)
variants <- variants[!(variants$clean == FALSE),]
variants <- variants[ , !(names(variants) %in% c("clean", "Flags"))]

# sort variants by position and get the position offset of each

gene.of.interest.row <- which(gene_file$name == gene.of.interest.symbol)
gene.of.interest.start <- as.numeric(strsplit(gene_file$start_position[gene.of.interest.row][[1]], split = "; ")[[1]][1])
gene.of.interest.end <- as.numeric(strsplit(gene_file$end_position[gene.of.interest.row][[1]], split = "; ")[[1]][1])
gene.of.interest.ch <- paste("chr",gene_file$chromosome[gene.of.interest.row],sep = "")
variants <- variants[order(variants$Position),]
# used to get the row number of the scores
get_position_offset <- function(n){
  offset <- variants$Position[n] - gene.of.interest.start + 1
  return(offset)
}
variants$distance.from.start <- sapply(1:dim(variants)[1], get_position_offset)

# find the cannonical exon (if any) each variant falls into
exon_regions <- strsplit(gene_file$exon[gene.of.interest.row][[1]], split = "; ")[[1]]
exon_regions <- lapply(exon_regions, function(x){substring(x, first = (as.integer(regexpr(":",x)[1])) + 1) })
exon_regions <- sapply(exon_regions, strsplit, split="-")
exon_regions <- sapply(exon_regions, as.numeric)
# correcting for different orientations
if(exon_regions[1,1] < exon_regions[2, length(exon_regions[1,])]){
  exon_regions <- exon_regions[,c(length(exon_regions[1,]):1)]
}

get_variant_exon <- function(position){
  closest_exon <- 0
  exon_dist <- abs(position - exon_regions[1,1])
  symb <- "+"
  for(i in 1:length(exon_regions[1,])){
    if((position >= exon_regions[1,i]) && (position <= exon_regions[2,i])){
      return(i)
    }
    if(abs(position - exon_regions[1,i]) < exon_dist){
      closest_exon <- i
      exon_dist <- abs(position - exon_regions[1,i])
      symb <- "+"
    }
    if(abs(position - exon_regions[2,i]) < exon_dist){
      closest_exon <- i
      exon_dist <- abs(position - exon_regions[2,i])
      symb <- "-"
    }
  }
  return(paste(closest_exon, symb, exon_dist, sep = ""))
}

variants$exon <- sapply(variants$Position,get_variant_exon)

# functions to get the ancestry of each variant in a nice format
get_ancestors <- function(n){
  ancestors = ""
  African <- variants$Allele.Count.African[n] - variants$Homozygote.Count.African[n]
  Ashkenazi_Jewish <- variants$Allele.Count.Ashkenazi.Jewish[n] - variants$Homozygote.Count.Ashkenazi.Jewish[n]
  East_Asian <- variants$Allele.Count.East.Asian[n] - variants$Homozygote.Count.East.Asian[n]
  European_Finnish <- variants$Allele.Count.European..Finnish.[n] - variants$Homozygote.Count.European..Finnish.[n]
  European <- variants$Allele.Count.European..Non.Finnish.[n] - variants$Homozygote.Count.European..Non.Finnish.[n]
  Latino <- variants$Allele.Count.Latino[n] - variants$Homozygote.Count.Latino[n]
  South_Asian <- variants$Allele.Count.South.Asian[n] - variants$Homozygote.Count.South.Asian[n]
  Other <- variants$Allele.Count.Other[n] - variants$Homozygote.Count.Other[n]
  ancestors <- paste(write_ancestors("African", African),write_ancestors("Ashkenazi Jewish", Ashkenazi_Jewish),
                     write_ancestors("East Asian", East_Asian), write_ancestors("Latino", Latino),
                     write_ancestors("European Finnish", European_Finnish),write_ancestors("European", European),
                     write_ancestors("South Asian", South_Asian),write_ancestors("Other", Other), sep = "")
  return(substr(ancestors, 1, nchar(ancestors) - 2))
}

write_ancestors <- function(name, number){
  ancestor <- ""
  if(length(number) > 0 && number != 0){
    ancestor <- paste(name," (", number, "); ", sep = "")
  }
  return(ancestor)
}

variants$ancestors <- sapply(1:dim(variants)[1], get_ancestors)
drops <- c("Allele.Count.African", "Allele.Number.African", "Homozygote.Count.African",
           "Allele.Count.Ashkenazi.Jewish", "Allele.Number.Ashkenazi.Jewish", "Homozygote.Count.Ashkenazi.Jewish",
           "Allele.Count.East.Asian", "Allele.Number.East.Asian", "Homozygote.Count.East.Asian",
           "Allele.Count.European..Finnish.", "Allele.Number.European..Finnish.", "Homozygote.Count.European..Finnish.",
           "Allele.Count.European..Non.Finnish.", "Allele.Number.European..Non.Finnish.", "Homozygote.Count.European..Non.Finnish.",
           "Allele.Count.Latino", "Allele.Number.Latino", "Homozygote.Count.Latino",
           "Allele.Count.South.Asian", "Allele.Number.South.Asian", "Homozygote.Count.South.Asian",
           "Allele.Count.Other", "Allele.Number.Other", "Homozygote.Count.Other")
variants <- variants[ , !(names(variants) %in% drops)]


# function to get gnomAD website for specific variant
get_gnomAD_website <- function(n){
  base_url <- "http://gnomad.broadinstitute.org/variant/"
  variantID <- paste(variants$Chrom[n], variants$Position[n], variants$Reference[n], variants$Alternate[n], sep = "-")
  return(paste(base_url,variantID,sep = ""))
}

variants$gnomAD.website <- sapply(1:dim(variants)[1], get_gnomAD_website)

#################################
# Get Genomic Scores
#################################

gr <- GRanges(seqnames=gene.of.interest.ch,
              IRanges(start=gene.of.interest.start:gene.of.interest.end, width=1))

# phastCons100way.UCSC.hg19 - phastCon scores are derived from the alignment of the human genome (hg19)
# and 99 other vertabrate species


gsco <- get0("gsco", ifnotfound = phastCons100way.UCSC.hg19)
citation(gsco) # the citation for the genomic scores
phastCon.scores <- scores(gsco, gr)

# used to get the phastCon score for each variant in the table
get_phastCon_score <- function(n){
  return(phastCon.scores[variants$distance.from.start[n]]$scores)
}
variants$phastCon.score <- sapply(1:dim(variants)[1], get_phastCon_score)

# fitCons.UCSC.hg19 - fitCons scores measure the fitness consequences of function annotation for the 
# human genome (hg19)

fitcon <- get0("fitcon", ifnotfound = fitCons.UCSC.hg19)
citation(fitcon) # the citation for the genomic scores
fitCon.scores <- scores(fitcon, gr)
# used to get the fitCon score for each variant in the table
get_fitCon_score <- function(n){
  return(fitCon.scores[variants$distance.from.start[n]]$scores)
}
variants$fitCon.score <- sapply(1:dim(variants)[1], get_fitCon_score)

# cadd.v1.3.hg19 - fitCons scores measure the fitness consequences of function annotation for the 
# human genome (hg19)
# These scores are rounded to provide faster lookup
if (!exists("cadd")){
  cadd <- getGScores("cadd.v1.3.hg19")
}
citation(cadd) # the citation for the genomic scores
# used to get the cadd score for each variant in the table
get_cadd_score <- function(n){
  if(variants$Annotation[n] != "missense"){
    return(NA)
  } 
  return(scores(cadd,GRanges(seqnames=gene.of.interest.ch,IRanges(
    start=variants$Position[n]:variants$Position[n], width=1)), ref=as.character(variants$Reference[n]),
    alt=as.character(variants$Alternate[n]))$scores)
}
variants$cadd.score <- sapply(1:dim(variants)[1], get_cadd_score)


# get mcap scores
if (!exists("mcap")){
  mcap <- getGScores("mcap.v1.0.hg19")
}
citation(mcap) # the citation for the genomic scores
# used to get the mcap score for each variant in the table
get_mcap_score <- function(n){
  if(variants$Annotation[n] != "missense"){
    return(NA)
  } 
  return(scores(mcap,GRanges(seqnames=gene.of.interest.ch,IRanges(
    start=variants$Position[n]:variants$Position[n], width=1)), ref=as.character(variants$Reference[n]),
    alt=as.character(variants$Alternate[n]))$scores)
}


variants$mcap.score <- sapply(1:dim(variants)[1], get_mcap_score)

variants$rest.api.url <- sapply(variants$Transcript.Consequence,function(x){
  paste("http://grch37.rest.ensembl.org/vep/human/hgvs/", gene.of.interest.symbol, ":",
        x ,"?content-type=application/json",sep="")})

variants$rest.api.info <- sapply(variants$rest.api.url, function(x){
  fromJSON(toJSON(content(GET(x))))})

variants$polyphen.score <- numeric(dim(variants)[1])
variants$polyphen.prediction <- character(dim(variants)[1])
variants$sift.score <- numeric(dim(variants)[1])
variants$sift.prediction <- character(dim(variants)[1])
variants$biotype <- character(dim(variants)[1])
variants$consequences.all <- character(dim(variants)[1])
variants$impact.all <- character(dim(variants)[1])

list_factors_pretty <- function(x){
  unique_factors <- unique(unlist(x))
  factor_string <- ""
  for(factor in unique_factors){
    factor_string <- paste(factor_string, factor, " (", sum(x == factor), "); ", sep = "")
  }
  return(factor_string)
}

extract_ensembl_api_info <- function(n){
  info <- variants$rest.api.info[[n]]$transcript_consequences[[1]]
  info <- info[info$gene_symbol == gene.of.interest.symbol,]
  variants$polyphen.score[n] <<- mean(as.numeric(unlist(info$polyphen_score)))
  variants$polyphen.prediction[n] <<- list_factors_pretty(info$polyphen_prediction)
  variants$sift.score[n] <<- mean(as.numeric(unlist(info$sift_score)))
  variants$sift.prediction[n] <<- list_factors_pretty(info$sift_prediction)
  variants$biotype[n] <<- list_factors_pretty(info$biotype)
  variants$consequences.all[n] <<- list_factors_pretty(info$consequence_terms)
  variants$impact.all[n] <<- list_factors_pretty(info$impact)
}

sapply(1:dim(variants)[1], extract_ensembl_api_info)

####################################################################
# You can change the score cutoffs below to match your research needs
####################################################################

# http://cadd.gs.washington.edu/info
# CADD minimum score:
score.cutoff.cadd <- 20

# https://link.springer.com/content/pdf/10.1186%2Fs40246-017-0104-8.pdf
# fitCons minimum score:
score.cutoff.fitCon <- 0.4

# phastCons minimum score:
score.cutoff.phastCon <- 0.55

# http://bejerano.stanford.edu/mcap/
# mcap minimum score:
score.cutoff.mcap <- 0.025 

# gnomAD maximum AF:
score.cutoff.AFgnomAD <- 0.01

# SIFT maximum score:
score.cutoff.sift <- 0.05

# PolyPhen2 mimimum score:
score.cutoff.polyphen <- 0.909

# functions to return the number of scores that pass/fail/NA

passes_cadd <- function(n){
  if(is.na(variants$cadd.score[n])){
    return(NA)
  }
  else if(variants$cadd.score[n] >= score.cutoff.cadd){
    return(1)
  }
  else{
    return(0)
  }
}

passes_fitCons <- function(n){
  if(is.na(variants$fitCon.score[n])){
    return(NA)
  }
  else if(variants$fitCon.score[n] >= score.cutoff.fitCon){
    return(1)
  }
}

passes_phastCons <- function(n){
  if(is.na(variants$phastCon.score[n])){
    return(NA)
  }
  else if(variants$phastCon.score[n] >= score.cutoff.phastCon){
    return(1)
  }
  else{
    return(0)
  }
}

passes_mcap <- function(n){
  if(is.na(variants$mcap.score[n])){
    return(NA)
  }
  else if(variants$mcap.score[n] >= score.cutoff.mcap){
    return(1)
  }
  else{
    return(0)
  }
}

passes_gnomAD <- function(n){
  if(is.na(variants$Allele.Frequency[n])){
    return(NA)
  }
  else if(variants$Allele.Frequency[n] <= score.cutoff.AFgnomAD){
    return(1)
  }
  else{
    return(0)
  }
}

passes_polyphen <- function(n){
  if(is.na(variants$polyphen.score[n])){
    return(NA)
  }
  else if(variants$polyphen.score[n] >= score.cutoff.polyphen){
    return(1)
  }
  else{
    return(0)
  }
}

passes_sift <- function(n){
  if(is.na(variants$sift.score[n])){
    return(NA)
  }
  else if(variants$sift.score[n] <= score.cutoff.sift){
    return(1)
  }
  else{
    return(0)
  }
}

overall_score <- function(n){
  pass <- 0
  fail <- 0
  n.a <- 0
  passed_list <- c(passes_gnomAD(n), passes_mcap(n), passes_phastCons(n), passes_fitCons(n), 
                   passes_cadd(n), passes_polyphen(n), passes_sift(n))
  for(i in passed_list){
    if(is.na(i)){
      n.a <- n.a + 1
    }
    else{
      pass <- pass + i
    }
  }
  fail <- length(passed_list) - pass - n.a
  return(c(pass, fail, n.a))
}

overall_score_lists <- lapply(1:dim(variants)[1], overall_score)
variants$num.pass <- sapply(1:dim(variants)[1], function(n){overall_score_lists[[n]][1]})
variants$num.fail <- sapply(1:dim(variants)[1], function(n){overall_score_lists[[n]][2]})
variants$num.na <- sapply(1:dim(variants)[1], function(n){overall_score_lists[[n]][3]})

# Erase unneeded rows
variants <- variants[ , !(names(variants) %in% c("rest.api.url", "rest.api.info", "source", "distance.from.start",
                                                 "rest.api.transcript.consequences","Consequence", "Filters...exomes", "Filters...genomes"))]

# sort variants by the number of scores that passed cutoffs
variants <- variants[order(-variants$num.pass),]

###################
# Write to CV
###################

write.csv(variants, file=paste(variants_file_path, variants_file_name, sep=""), row.names = FALSE)

##################
# Make Plot
##################

ensembl <- get0("ensembl", ifnotfound = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)) 
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome =  gene.of.interest.ch)
gene.image <- makeGene(id = gene.of.interest.symbol, type = "hgnc_symbol", biomart = ensembl) 
genomeAxis <- makeGenomeAxis(add53 = TRUE, add35=TRUE)
expres <- makeGenericArray(intensity = as.matrix(variants$num.pass), probeStart = as.numeric(
  variants$Position), dp = DisplayPars(type = "dot", lwd = 2, pch = "o")) # shows number of scores passed
gdPlot(list(genomeAxis, Exons=gene.image, "Scores Passed"= expres),
       minBase = gene.of.interest.start, maxBase =gene.of.interest.end, labelCex = 2) # plots all 3 images
