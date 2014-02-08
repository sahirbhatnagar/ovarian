source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)
edgeRUsersGuide()

setkey(final.data)

#both of these work
dt <- final.data

#dt[, genesum:=apply(dt[,-1, with=FALSE],1, sum)]
#this one is faster
dt[, list(sum=rowSums(dt[, -1, with=FALSE]),
          mean=rowMeans(dt[, -1, with=FALSE])]
head(dt)
dim(dt)
head(final.data)

boxplot(dt[, -1, with=FALSE])
dt[J(1)]

example(data.table)
