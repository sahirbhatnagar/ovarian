##################################
# R source code file used to find people who have both phenotype and RNA_seq data
# Created by Sahir Rai Bhatnagar, Jan 17, 2014
# Updated Jan 17, 2014
# hosted on Github repo 'sahirbhatnagar/ovarian'
##################################

# Find people for who we have both phenotype and gene expression data --------------------
setwd(gene.path)
file.names <- list.files(pattern = "[.]txt$")
k <- lapply(o.pheno2$ID, function(i) grep(i,file.names, perl=TRUE, value=TRUE))



# files for which we have phenotype data and gene expression data (n=133)
matches <- unlist(k)

setwd(phenotype.path)
