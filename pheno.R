##################################
# R source code file used to remove missing data from phenotype file 'TCGA_OV_PhenotypeFile.csv'
# Created by Sahir Rai Bhatnagar
# Updated Jan 8, 2014
# hosted on Github repo 'sahirbhatnagar/ovarian'
##################################


# Import Data -------------------------------------------------------------
getwd()
setwd("/home/sahir/git_repositories/ovarian/")

# From Github
# library(repmis)
# url<-paste0("https://raw.github.com/sahirbhatnagar/","ovarian/master/","TCGA_OV_PhenotypeFile.csv")
# o.pheno<-source_data(url)

# or locally
o.pheno <- read.csv("TCGA_OV_PhenotypeFile.csv", header = T)[, c("ID", "Sample", 
    "Barcode.SNP", "TP53.class", "AgeAtDiagnosis..yrs.", "TUMORSTAGE", "TUMORGRADE", 
    "ProgressionFreeStatus", "ProgressionFreeSurvival..mos..")]


# Replace spaces with NA's ------------------------------------------------
o.pheno[o.pheno == ""] <- NA
o.pheno[o.pheno == "Missing"] <- NA
o.pheno[o.pheno == "TP53 missense (2)"] <- "TP53 missense"
o.pheno[o.pheno == "TP53 null/missense"] <- "TP53 null"


# Delete missing values ---------------------------------------------------
o.pheno2 <- o.pheno[complete.cases(o.pheno), ]
o.pheno2$TP53.class <- factor(o.pheno2$TP53.class)
o.pheno2$TUMORSTAGE <- factor(o.pheno2$TUMORSTAGE)
o.pheno2$TUMORGRADE <- factor(o.pheno2$TUMORGRADE)


# Create Survival Variables -----------------------------------------------
o.pheno2$cens[o.pheno2$ProgressionFreeStatus == "Recurred/Progressed"] <- 1
o.pheno2$cens[o.pheno2$ProgressionFreeStatus == "DiseaseFree"] <- 0
o.pheno2$time <- as.numeric(as.character(o.pheno2$ProgressionFreeSurvival..mos..))


# Convert tumorstage into numerical ----------------------------------------
o.pheno2$tumorstage <- sapply(o.pheno2$TUMORSTAGE, function(x) {
    if (x %in% c("IIA", "IIB", "IIC")) 
        2 else if (x %in% c("IIIA", "IIIB", "IIIC")) 
        3 else 4
}) 


# Find people for who we have both phenotype and gene expression data --------------------
setwd("/home/sahir/git_repositories/ovarian/data")
file.names <- list.files(pattern = "[.]txt$")
k <- lapply(o.pheno2$ID, function(i) grep(i,file.names, perl=TRUE, value=TRUE))

# files for which we have phenotype data and gene expression data (n=133)
matches <- unlist(k)


# Create list of gene names -----------------------------------------------
#import gene expression file
gene.exp <- read.table(matches[1], header=TRUE)

#extract gene names
gene.exp$gene <- gsub("\\|", " ", as.character(gene.exp$gene))
gene.exp$gene <- sub(" .*", "", gene.exp$gene)

#remove unknown genes
gene.exp<-gene.exp[!(gene.exp$gene=="?"),]

#unique gene names
unique.genes <- unique(gene.exp$gene)



# Create function to extract gene names for each file ---------------------

ex.genes <- function(filename, type="raw_counts"){
  #import gene expression file
  gene.exp <- read.table(filename, header=TRUE)
  
  #extract gene names
  gene.exp$gene <- gsub("\\|", " ", as.character(gene.exp$gene))
  gene.exp$gene <- sub(" .*", "", gene.exp$gene)
  
  #remove unknown genes
  gene.exp<-gene.exp[!(gene.exp$gene=="?"),]
  
  #unique gene names
  unique.genes <- unique(gene.exp$gene)
  
  return(gene.exp)
}

d1<-ex.genes(matches[1])
d2<-ex.genes(matches[2])

k<-sum(!(d1$gene == unique.genes))


a<-c(1,2,3,4,5)
b<-c()

# Merge all gene expresssion files into one text file ---------------------

library(doParallel)
registerDoParallel(cores = 4)
library(foreach)
all.gene <- foreach(i = matches, .combine = rbind) %dopar% read.table(i, header = TRUE)

#rm(gene.exp)
rm(all.gene)
