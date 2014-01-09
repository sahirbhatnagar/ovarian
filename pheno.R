##################################
# R source code file used to remove missing data from phenotype file 'TCGA_OV_PhenotypeFile.csv'
# Created by Sahir Rai Bhatnagar
# Updated Jan 8, 2014
# hosted on Github repo 'sahirbhatnagar/ovarian'
##################################


# Import Data -------------------------------------------------------------
getwd()
setwd("/home/sahir/git_repositories/ovarian")

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


# Conver tumorstage into numerical ----------------------------------------
o.pheno2$tumorstage <- sapply(o.pheno2$TUMORSTAGE, function(x) {
    if (x %in% c("IIA", "IIB", "IIC")) 
        2 else if (x %in% c("IIIA", "IIIB", "IIIC")) 
        3 else 4
}) 


charmatch
