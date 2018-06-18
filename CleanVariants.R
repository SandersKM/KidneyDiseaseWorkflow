
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
variants <- variants[ , !(names(variants) == "clean")]

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



