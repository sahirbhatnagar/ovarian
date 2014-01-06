############################################################################
# expression data
library(plyr)
getwd()
setwd(paste(getwd(),"/data", sep=""))

file.names<-list.files(pattern = "[.]txt$")
names(file.names)<-basename(file.names)

#read in files
rnaseq<-ldply(file.names[1], read.table, header=TRUE)[,-1]

#extract gene names
rnaseq$pp<-gsub("\\|"," ", as.character(rnaseq$gene))
rnaseq$pp<-sub(" .*", "", rnaseq$pp)
rnaseq$pp<-factor(rnaseq$pp)

str(rnaseq)

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



######
#foreach function for loops on multiple cores
library(foreach)
library(doMC)

registerDoMC(cores=4)

getDoParWorkers()
getDoParName()
getDoParVersion()





#  x <- iris[which(iris[,5] != "setosa"), c(1,5)]
#  trials <- 10000
#  ptime <- system.time({
#    r <- foreach(icount(trials), .combine=cbind) %dopar% {
#      ind <- sample(100, 100, replace=TRUE)
#      result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
#     coefficients(result1)
#      }
#    })[3]
#  ptime

library(plyr)
library(bigmemory)
library(biganalytics)
getwd()
setwd(paste(getwd(),"/data", sep=""))

file.names<-list.files()
names(file.names)<-basename(file.names)

#read in files using read.big.matrix
#cannot have character type in read.big.matrix... type="char"
# The 'char' type is a C++ data type used to store integer 
# values that represent ASCII character codes (a single character, not a character string).
# To store character strings in a big.matrix, you'll have to re-code the strings as numeric 
# values (or convert to factors, then to numeric values).
# If you need to store character data in a very large data set, you may want to look 
# into the 'ff' package. In my experience it has a steep learning curve and the 
# documentation is somewhat lacking, but it does have that functionality.
# 
# For further details on dealing with large data sets, you can check out the CRAN Task 
# View here: http://cran.r-project.org/web/views/HighPerformanceComputing.html

x<-read.big.matrix(file.names[1], header=TRUE, backingfile="x.bin", 
                   sep="\t", descriptorfile="x.desc", type="char")

summary(x)
x[155:160,]
head(x)

#extract gene names
x$pp<-gsub("\\|"," ", as.character(x$gene))
x$pp<-sub(" .*", "", x$pp)




# combine all files
files <- list.files(pattern = "[.]txt$")[1:5]

library(rbenchmark)
benchmark(replications = 5, order = "user.self",
          LAPPLY = {
            read.all <- lapply(files, read.table, header = TRUE)
            data1 <- do.call(rbind, read.all)
          },
          FOREACH = {
            library(doParallel)
            registerDoParallel(cores = 4)
            library(foreach)
            data2 <- foreach(i = files, .combine = rbind) %dopar% read.table(i, header = TRUE)
          }
)

library(compare)
all.equal(data1, data2)




