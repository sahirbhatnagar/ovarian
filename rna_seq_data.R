############################################################################
# expression data
library(plyr)
getwd()

file.names<-list.files("data")
names(file.names)<-basename(file.names)

#read in files
rnaseq<-ldply(file.names[1], read.table, header=TRUE)[,-1]

#extract gene names
rnaseq$pp<-gsub("\\|"," ", as.character(rnaseq$gene))
rnaseq$pp<-sub(" .*", "", rnaseq$pp)


########
#sqldf

library(sqldf)
sqldf("select * from rnaseq")
#hello world 


library ( devtools )
#install_github(repo = "pbdMPI" , username = "RBigData")
#install_github(repo = "pbdSLAP" , username = "RBigData")
install_github(repo = "pbdNCDF4" , username = "RBigData")
install_github(repo = "pbdBASE" , username = "RBigData")
#install_github(repo = "pbdDMAT", username = "RBigData")
install_github(repo = "pbdDEMO" , username = "RBigData")


#####
# bigmemory
library(bigmemory)
library(biganalytics)

