##################################
# R source code file used to create data frame of RNA_seq data
# Created by Sahir Rai Bhatnagar, Jan 17, 2014
# Updated Jan 17, 2014
# hosted on Github repo 'sahirbhatnagar/ovarian'
##################################

getwd()
setwd("/home/sahir/git_repositories/ovarian/")
source("pheno.R")
source("matches.R")



# Create gene expression data frame for given id number -----------------------------------------------
# import gene expression file
gene.exp <- read.table(matches[10], header=TRUE)

#extract gene names
gene.exp$gene <- gsub("\\|", " ", as.character(gene.exp$gene))
gene.exp$gene <- sub(" .*", "", gene.exp$gene)

#remove unknown genes
gene.exp<-gene.exp[!(gene.exp$gene=="?"),]

#rename raw_counts column to ID number of person
type <- "raw_counts"
id.number.original <- unique(substring(gene.exp$barcode,first=1,last=12))

id.number<-gsub("-",c("_"),id.number.original)

colnames(gene.exp)[grepl(type,colnames(gene.exp))] <- id.number

#remove all columns except expression counts and gene names
gene.exp <- gene.exp[,c("gene",id.number)]

# Now we need to average the measurements across each gene ----------------
library(data.table)
gene.exp.big <- as.data.table(gene.exp)
gene.exp.avg <- gene.exp.big[ , lapply(.SD, mean),by='gene', .SDcols=id.number ]

# Create pheno data frame for given id number -----------------------------

#bring in phenotype data for a given ID number
j <- o.pheno2[o.pheno2$ID==id.number.original,c("ID","TP53.class","AgeAtDiagnosis..yrs.","cens","time","tumorstage")]

library(reshape)
pheno.final <- melt(j, id=c("ID" ), measure.vars=c("AgeAtDiagnosis..yrs.","cens","time","tumorstage"))[,2:3]
colnames(pheno.final)<-c("gene",id.number)


# Merge gene exp with pheno data frames -----------------------------------
k <- rbind(gene.exp,pheno.final)


# Create function to merge gene exp with pheno data ---------------------

library(data.table)
library(reshape)
f.gene <- function(filename, type="raw_counts"){
  #import gene expression file
  gene.exp <- read.table(filename, header=TRUE)
  
  #extract gene names
  gene.exp$gene <- gsub("\\|", " ", as.character(gene.exp$gene))
  gene.exp$gene <- sub(" .*", "", gene.exp$gene)
  
  #remove unknown genes
  gene.exp<-gene.exp[!(gene.exp$gene=="?"),]
  
  #rename raw_counts column to ID number of person
  id.number.original <- unique(substring(gene.exp$barcode,first=1,last=12))
  
  id.number <- gsub("-",c("_"),id.number.original)
  
  colnames(gene.exp)[grepl(type,colnames(gene.exp))] <- id.number
  
  #remove all columns except expression counts and gene names
  gene.exp <- gene.exp[,c("gene",id.number)]
  
  # Now we need to average the measurements across each gene ----------------
  gene.exp <- as.data.table(gene.exp)
  gene.exp <- gene.exp[ , lapply(.SD, mean),by='gene', .SDcols=id.number ]
  
  #bring in phenotype data for a given ID number
  j <- o.pheno2[o.pheno2$ID==id.number.original,c("ID","TP53.class","AgeAtDiagnosis..yrs.","cens","time","tumorstage")]
  pheno.final <- melt(j, id=c("ID" ), measure.vars=c("AgeAtDiagnosis..yrs.","cens","time","tumorstage"))[,2:3]
  colnames(pheno.final)<-c("gene",id.number)
  
  # Merge gene exp with pheno data frames -----------------------------------
  k<-rbind(gene.exp,pheno.final)
  
  return(k)
}


# merge two data.tables. this works
l <- f.gene(matches[2], type="raw_counts")
m <- f.gene(matches[3], type="raw_counts")

setkey(l, gene)
res <- l[m]








# Merge all gene expresssion files into one text file ---------------------

library(doParallel)
registerDoParallel(cores = 4)
library(foreach)
all.gene <- foreach(i = matches, .combine = rbind) %dopar% read.table(i, header = TRUE)

#rm(gene.exp)
rm(all.gene)
