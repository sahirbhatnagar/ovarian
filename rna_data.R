##################################
# R source code file of the function 'f.gene' which created data frame of RNA_seq and pehno
# Created by Sahir Rai Bhatnagar, Jan 17, 2014
# Updated Jan 17, 2014
# hosted on Github repo 'sahirbhatnagar/ovarian'
##################################

# Create function to merge gene exp with pheno data ---------------------

library(data.table)
library(reshape)
f.gene <- function(filename, type="raw_counts", dir="/home/sahir/git_repositories/ovarian/data"){
  setwd(dir)
  #import gene expression file
  gene.exp <- read.table(filename, header=TRUE)
  
  #extract gene names
  gene.exp$gene <- gsub("\\|", " ", as.character(gene.exp$gene))
  gene.exp$gene <- sub(" .*", "", gene.exp$gene)
  
  #remove unknown genes
  gene.exp<-gene.exp[!(gene.exp$gene=="?"),]
  
  #rename raw_counts column to ID number of person
  id.number.original <- unique(substring(gene.exp$barcode,first=1,last=16))
  
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


# Merge all gene expresssion files into one data file ---------------------

library(doParallel)
registerDoParallel(cores = 4)
library(foreach)
#system.time(res <- lapply(matches, f.gene))
#this is faster than lapply
all.gene <- foreach(i = matches) %dopar% f.gene(i)



# This produces the final dataset to be analysed --------------------------

lapply(all.gene, function(i) setkey(i,gene))
final.data <- Reduce(merge,all.gene)










