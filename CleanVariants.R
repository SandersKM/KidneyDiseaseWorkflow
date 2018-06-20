suppressPackageStartupMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(GenomicScores))
suppressPackageStartupMessages(library(phastCons100way.UCSC.hg19))
suppressPackageStartupMessages(library(fitCons.UCSC.hg19))
suppressPackageStartupMessages(library(cadd.v1.3.hg19))
suppressPackageStartupMessages(library(mcap.v1.0.hg19))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library('GenomicFeatures'))


samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                          package="GenomicFeatures")
txdb <- loadDb(samplefile)
gene.of.interest.ch <- paste("chr",gene_file$chromosome[gene.of.interest.row],sep = "")
seqlevels(txdb) <- gene.of.interest.ch
# to reset seqlevel use seqlevels(txdb) <- seqlevels0(txdb)
#TODO - figure out the exon thing :-/

# Change the symbol to correspond to the gene of interest
gene.of.interest.symbol <- "MUC1"
gene.of.interest.row <- which(gene_file$Symbol == "MUC1")

# download the variant CSV from the gnomAD browser for the 
# gene of interest. Write the path file below
all_variants_path <- "/Users/ksanders/Downloads/MUC1_genes_all.csv"
variants <- read.csv(all_variants_path, header = TRUE, sep=",")

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

# TODO - label exon numbers (maybe get exons from somewhere?)

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
  if(number != 0){
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

# used to get the row number of the scores
get_position_offset <- function(n){
  offset <- variants$Position[n] - gene.of.interest.start + 1
  return(offset)
}

gene.of.interest.start <- gene_file$start_position_on_the_genomic_accession[gene.of.interest.row]
gene.of.interest.end <- gene_file$end_position_on_the_genomic_accession[gene.of.interest.row]
variants$distance.from.start <- sapply(1:dim(variants)[1], get_position_offset)

# phastCons100way.UCSC.hg19 - phastCon scores are derived from the alignment of the human genome (hg19)
# and 99 other vertabrate species

gsco <- phastCons100way.UCSC.hg19
citation(gsco) # the citation for the genomic scores
phastCon.scores <- scores(gsco, gr)

# used to get the phastCon score for each variant in the table
get_phastCon_score <- function(n){
  return(phastCon.scores[variants$distance.from.start[n]]$scores)
}
variants$phastCon.score <- sapply(1:dim(variants)[1], get_phastCon_score)

# fitCons.UCSC.hg19 - fitCons scores measure the fitness consequences of function annotation for the 
# human genome (hg19)

fitcon <- fitCons.UCSC.hg19
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
cadd <- getGScores("cadd.v1.3.hg19")
citation(cadd) # the citation for the genomic scores
# used to get the cadd score for each variant in the table
get_cadd_score <- function(n){
  if(variants$Annotation[n] != "missense"){
    return("NA")
  } 
  return(scores(cadd,GRanges(seqnames=gene.of.interest.ch,IRanges(
    start=variants$Position[n]:variants$Position[n], width=1)), ref=as.character(variants$Reference[n]),
    alt=as.character(variants$Alternate[n]))$scores)
}
variants$cadd.score <- sapply(1:dim(variants)[1], get_cadd_score)


# ah = AnnotationHub()
# ah <- subset(ah, species == "Homo sapiens")
# ah
# grch37 <- query(ah, "GRCh37", "hg19")
# grch37$genome
# display(grch37)




# get mcap scores
mcap <- getGScores("mcap.v1.0.hg19")
citation(mcap) # the citation for the genomic scores
# used to get the mcap score for each variant in the table
get_mcap_score <- function(n){
  if(variants$Annotation[n] != "missense"){
    return("NA")
  } 
  return(scores(mcap,GRanges(seqnames=gene.of.interest.ch,IRanges(
    start=variants$Position[n]:variants$Position[n], width=1)), ref=as.character(variants$Reference[n]),
    alt=as.character(variants$Alternate[n]))$scores)
}
variants$mcap.score <- sapply(1:dim(variants)[1], get_mcap_score)


