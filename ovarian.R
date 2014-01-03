library(survival)
setwd("/home/sahir/Dropbox/PhD/Year 1/Prelim Analysis Nov 2013/ovarian/")

o.pheno<-read.csv("/home/sahir/Dropbox/PhD/Year 1/Prelim Analysis Nov 2013/ovarian/TCGA_OV_PhenotypeFile.csv",
                  header=T)[,c("ID","Sample","Barcode.SNP","TP53.class","AgeAtDiagnosis..yrs.",
                                                         "TUMORSTAGE","TUMORGRADE","ProgressionFreeStatus",
                                                         "ProgressionFreeSurvival..mos..")]
o.pheno[o.pheno==""]<-NA
o.pheno[o.pheno=="Missing"]<-NA
o.pheno[o.pheno=="TP53 missense (2)"]<-"TP53 missense"
o.pheno[o.pheno=="TP53 null/missense"]<-"TP53 null"

#o.pheno2<-o.pheno[complete.cases(o.pheno[,c("TP53.class","ProgressionFreeStatus","ProgressionFreeSurvival..mos..")]),  ]
o.pheno2<-o.pheno[complete.cases(o.pheno),  ]
o.pheno2$cens[o.pheno2$ProgressionFreeStatus=="Recurred/Progressed"]<-1
o.pheno2$cens[o.pheno2$ProgressionFreeStatus=="DiseaseFree"]<-0
o.pheno2$time<-as.numeric(as.character(o.pheno2$ProgressionFreeSurvival..mos..))
o.pheno2$TP53.class<-factor(o.pheno2$TP53.class)
o.pheno2$TUMORSTAGE<-factor(o.pheno2$TUMORSTAGE)
o.pheno2$TUMORGRADE<-factor(o.pheno2$TUMORGRADE)


#conver tumorstage into numerical
o.pheno2$tumorstage<-sapply(o.pheno2$TUMORSTAGE, function(x) {
  if (x %in% c("IIA","IIB","IIC")) 2 else if 
     (x %in% c("IIIA","IIIB","IIIC")) 3 else 4 } )

table(o.pheno2$tumorstage)
table(o.pheno2$TUMORSTAGE)
str(o.pheno2)

#create survival object
pheno.obj<-with(o.pheno2,  Surv(time,cens))

kmfit<-with(o.pheno2, survfit(pheno.obj~TP53.class))
summary(kmfit)
source("/home/sahir/Dropbox/PhD/Year 1/Prelim Analysis Nov 2013/ovarian/ggsurv.R")
pl2<-ggsurv(kmfit, xlab='Time (months)')
med.surv <- data.frame(time = c(16.6,16.6, 16.8,16.8), quant = c(.5,0,.5,0),
                       TP53.class = c('TP53 missense', 'TP53 missense', 'TP53 null', 'TP53 null'))
pl2 + geom_line(data = med.surv, aes(time, quant, group = TP53.class), 
                col = 'darkblue', linetype = 3) +
  geom_point(data = med.surv, aes(time, quant, group =TP53.class), col = 'darkblue')+
  annotate("text", x = 14, y = 0.25, label = "median")

#log-rank test
survdiff(pheno.obj~TP53.class, data=o.pheno2)

#wilcoxon test
wilcox.test(time~TP53.class, data=o.pheno2)

#cox-model
library(rms)
fit.cox<-cph(pheno.obj ~ AgeAtDiagnosis..yrs.+TUMORGRADE+tumorstage+strat(TP53.class), 
             data=o.pheno2, x=TRUE, y=TRUE)

describe(o.pheno2[,c("AgeAtDiagnosis..yrs.","TUMORGRADE","tumorstage")])


#parametric model
with(o.pheno2, survplot(survfit(pheno.obj ~TP53.class), conf='none', fun=qnorm, logt=TRUE  ))

fit.aft<-psm(pheno.obj~AgeAtDiagnosis..yrs.+TUMORGRADE+tumorstage+ TP53.class, data=o.pheno2, dist='lognormal')

r<-resid(fit.aft)

with(o.pheno2, survplot(r, TP53.class, label.curve=FALSE))
with(o.pheno2, survplot(r, AgeAtDiagnosis..yrs., label.curve=FALSE))
with(o.pheno2, survplot(r, TUMORGRADE, label.curve=FALSE))
with(o.pheno2, survplot(r, tumorstage, label.curve=FALSE))

Function(fit.cox)

latex(fit.cox, file="eqn.tex")

latex(print(fit.cox))

plot(fit.cox)
     
gendata(fit.aft, AgeAtDiagnosis..yrs.=50, TUMORGRADE="G3", tumorstage=4)

fit.aft<-psm(pheno.obj~AgeAtDiagnosis..yrs.+TUMORGRADE+tumorstage, data=o.pheno2 )

survplot(fit.aft, AgeAtDiagnosis..yrs.=50, TUMORGRADE="G3", tumorstage=4 )


############################################################################
# expression data
library(plyr)
getwd()
setwd("/home/sahir/Dropbox/PhD/Year 1/Prelim Analysis Nov 2013/ovarian/RNASeq/")
file.names<-list.files()
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