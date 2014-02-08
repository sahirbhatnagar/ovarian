##################################
# R source makefile for RNA_seq ovarian data which produces the final data set to be analysed
# Created by Sahir Rai Bhatnagar, Jan 17, 2014
# Updated Jan 17, 2014
# hosted on Github repo 'sahirbhatnagar/ovarian'
##################################

# path for the phenotype file
phenotype.path <- "/home/sahir/git_repositories/ovarian/"
#path for the gene expression data
gene.path <- "/home/sahir/git_repositories/ovarian/data"
type <- "RPKM"
# type <- "median_length_normalized"
# type <- "raw_counts"

source("pheno.R")
source("matches.R")
source("rna_data.R")

getwd()
write.csv(final.data, "finaldata")
f<-read.csv("finaldata", header=TRUE)
