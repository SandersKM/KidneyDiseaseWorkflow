
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

