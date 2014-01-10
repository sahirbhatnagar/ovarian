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


# Create gene expression data frame for given id number -----------------------------------------------
#import gene expression file
gene.exp <- read.table(matches[10], header=TRUE)

#extract gene names
gene.exp$gene <- gsub("\\|", " ", as.character(gene.exp$gene))
gene.exp$gene <- sub(" .*", "", gene.exp$gene)

#remove unknown genes
gene.exp<-gene.exp[!(gene.exp$gene=="?"),]

#rename raw_counts column to ID number of person
type <- "raw_counts"
id.number <- unique(substring(gene.exp$barcode,first=1,last=12))

id.number<-gsub("-",c("_"),id.number)

colnames(gene.exp)[grepl(type,colnames(gene.exp))] <- id.number

#remove all columns except expression counts and gene names
gene.exp <- gene.exp[,c("gene",id.number)]



library(data.table)
gene.exp.big <- as.data.table(gene.exp)
str(gene.exp.big)

gene.exp.avg <- gene.exp.big[,mean(id.number),by='gene']


# Create pheno data frame for given id number -----------------------------

#bring in phenotype data for a given ID number
j<-o.pheno2[o.pheno2$ID==id.number,c("ID","TP53.class","AgeAtDiagnosis..yrs.","cens","time","tumorstage")]

library(reshape)
pheno.final <- melt(j, id=c("ID" ), measure.vars=c("AgeAtDiagnosis..yrs.","cens","time","tumorstage"))[,2:3]
colnames(pheno.final)<-c("gene",id.number)






# Merge gene exp with pheno data frames -----------------------------------
k<-rbind(gene.exp,pheno.final)


# Create function to merge gene exp with pheno data ---------------------

f.gene <- function(filename, type="raw_counts"){
  #import gene expression file
  gene.exp <- read.table(filename, header=TRUE)
  
  #extract gene names
  gene.exp$gene <- gsub("\\|", " ", as.character(gene.exp$gene))
  gene.exp$gene <- sub(" .*", "", gene.exp$gene)
  
  #remove unknown genes
  gene.exp<-gene.exp[!(gene.exp$gene=="?"),]
  
  #rename raw_counts column to ID number of person
  id.number <- unique(substring(gene.exp$barcode,first=1,last=12))
  colnames(gene.exp)[grepl(type,colnames(gene.exp))] <- id.number
  
  #remove all columns except expression counts and gene names
  gene.exp <- gene.exp[,c("gene",id.number)]
  
  gene.exp.avg <- ddply(gene.exp,.(gene), function(x) data.frame(gene=x$gene[1],
                                                                 id.number=mean(x[,2])))
  
  #bring in phenotype data for a given ID number
  j<-o.pheno2[o.pheno2$ID==id.number,c("ID","TP53.class","AgeAtDiagnosis..yrs.","cens","time","tumorstage")]
  pheno.final <- melt(j, id=c("ID" ), measure.vars=c("AgeAtDiagnosis..yrs.","cens","time","tumorstage"))[,2:3]
  colnames(pheno.final)<-c("gene",id.number)
  
  # Merge gene exp with pheno data frames -----------------------------------
  k<-rbind(gene.exp.avg,pheno.final)
    
  return(k)
}


x<-f.gene(matches[2])[1:50,]
y<-f.gene(matches[3])[1:50,]

dat<-merge(x,y, by.x="gene", all=T)



# Merge all gene expresssion files into one text file ---------------------

library(doParallel)
registerDoParallel(cores = 4)
library(foreach)
all.gene <- foreach(i = matches, .combine = rbind) %dopar% read.table(i, header = TRUE)

#rm(gene.exp)
rm(all.gene)
