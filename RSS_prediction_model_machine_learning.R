# Author :  Jun Chen
# date :  Aug,8 2013
# Goal :  convert data format
# 
# Sex check theory.
# A PROBLEM arises if the two sexes do not match, or if the SNP data or pedigree data 
# are ambiguous with regard to sex. A male call is made if F is more than 0.8; a femle call is made if F is less than 0.2. 
# F is The actual Z chromosome inbreeding (homozygosity) estimate


# Method 1.  
# Method 2.   
# Method 3.   
# Method 4.  
# Method 5.   
# Method 6.   



# setwd("/home/chenju/project_RSS/RSS") # change here
#setwd("X:/project_RSS/RSS") # change here
setwd("Z:/project_RSS/RSS") # change here

#snp_60k <- read.table("/work/d/Jun/project_gs/imputation/data/reference/SNP_60K_Sorted.txt",header=TRUE,sep="\t")
#snp_60k <- read.table("W:/Jun/project_gs/imputation/data/reference/SNP_60K_Sorted.txt",header=TRUE,sep="\t")
snp_60k <- read.table("V:/Jun/project_gs/imputation/data/reference/SNP_60K_Sorted.txt",header=TRUE,sep="\t")
snp_60k_chrZ <- subset(snp_60k, Chr=="Z") # This is for sex checking.
snp_60k_chrW <- subset(snp_60k, Chr=="W")

# file1
#temp <- read.table(file="/home/chenju/project_RSS/RSS/12-163_20120928_190 samples_sim-rss_FinalReport_Matrix AB.txt",head=TRUE,sep="\t",check.names=F,colClasses="character",na.strings = c("--"),skip=9, row.names=1)
#temp <- read.table(file="X:/project_RSS/RSS/12-163_20120928_190 samples_sim-rss_FinalReport_Matrix AB.txt",head=TRUE,sep="\t",check.names=F,colClasses="character",na.strings = c("--"),skip=9, row.names=1)
temp <- read.table(file="Z:/project_RSS/RSS/12-163_20120928_190 samples_sim-rss_FinalReport_Matrix AB.txt",head=TRUE,sep="\t",check.names=F,colClasses="character",na.strings = c("--"),skip=9, row.names=1)
Line_SNP <- t(temp)
num_4 <- match(as.character(unlist(snp_60k[,1])),colnames(Line_SNP))
num_5 <- which(num_4>0)
length(num_5)
Mydata_temp <- Line_SNP[,c(num_4)]


# file 2 for sex checking.
num_4 <- match(as.character(unlist(snp_60k_chrZ[,1])),colnames(Line_SNP))
num_5 <- which(num_4>0)
length(num_5)
Mydata_temp_sex <- Line_SNP[,c(num_4)]

Sex_results <- array(0,dim(Mydata_temp_sex)[1])
Total <- dim(Mydata_temp_sex)[2]
for(i in 1:dim(Mydata_temp_sex)[1])
{
Hom <- length(which(Mydata_temp_sex[i,] == "AA" | Mydata_temp_sex[i,] == "BB"))
Na_count <-length(which(as.character(Mydata_temp_sex[i,])== "NA"))
Hom_rate <- Hom/(Total-Na_count)
if(Hom_rate >= 0.95) Sex_results[i] <- "F"  # Het -> Male (ZZ),  Hom -> Female (ZW) 
if(Hom_rate < 0.95 )  Sex_results[i] <- "M"   
}

rownames(Sex_results) <- rownames(Mydata_temp_sex)
write.table(Sex_results,file="Sex_Est_Jun.txt",append=FALSE,quote=FALSE,sep=" ",row.names=TRUE,col.names=TRUE)
# end of sex estimation.

# Summary of SNP
Mydata_temp_status <- substr(rownames(Mydata_temp),1,3)
Phe <- Mydata_temp_status
myData <- cbind(Phe,Mydata_temp)
library(SNPassoc)
DatSNP<-setupSNP(data=myData,colSNPs=2:dim(myData)[2],sep="")
DatSNP$Phe <- as.factor(DatSNP$Phe)
results <- summary(DatSNP)   # good summary program.
write.table(results,file="RSS_SNP_chips_summary.txt",append=FALSE,quote=FALSE,sep=" ",row.names=TRUE,col.names=TRUE)

# model = c("codominant", "dominant", recessive, overdominant, log-additive or all),quantitative = FALSE
Association_results <- WGassociation(Phe, data=DatSNP, model ="recessive",level = 0.95)

# Pick up top 50, and classification and clust.
# End of here.

data(SNPs.info.pos)  # marker file
> myData.o<-setupSNP(SNPs, colSNPs=6:40, sort=TRUE,
+ info=SNPs.info.pos, sep="")
association(casco~snp10001, data=myData, model=c("cod","log"))

# **********************************************************************************
# Machine learning from here.


#source("http://bioconductor.org/biocLite.R")
#biocLite("genefilter")
#biocLite("Biobase")
#biocLite("bioDist")
#biocLite("MLInterfaces")
#install.packages("RColorBrewer")
#install.packages("randomForest")

library("Biobase")
library("genefilter")
library("bioDist")
library("MLInterfaces")


neg = which(Mydata_temp_status == "neg")
pos = which(Mydata_temp_status == "pos")
S1 = sample(neg, 50, replace=FALSE)
S2 = sample(pos, 50, replace = FALSE)
TrainInd = c(S1, S2)
TestInd = setdiff(1:190, TrainInd)
#f2 <- factor(Mydata_temp_status)



Mydata <- Mydata_temp

genotype11<-"AA"
genotype12<-"AB"
genotype22<-"BB"
for(i in 1:dim(Mydata)[2])
{
  if(length(which(Mydata[,i]==as.character(genotype11)))>0) Mydata[c(which(Mydata[,i]==as.character(genotype11))),i]="0"
  if(length(which(Mydata[,i]==as.character(genotype12)))>0) Mydata[c(which(Mydata[,i]==as.character(genotype12))),i]="1"
  if(length(which(Mydata[,i]==as.character(genotype22)))>0) Mydata[c(which(Mydata[,i]==as.character(genotype22))),i]="2"
}


for(i in 1:dim(Mydata)[2])
{
  Mydata[,i] <- as.numeric(Mydata[,i])
}

Phe <- Mydata_temp_status

# cluster here, see what happened.






# Traintt = rowttests(Mydata[,TrainInd],f2[TrainInd])
#a <- table(Phe[TrainInd],Mydata[i,TrainInd])
#ordTT = order(abs(Traintt$statistic), decreasing=TRUE)


Mydata2 <- cbind(Phe,Mydata)
fNtt = Mydata2[,c(1,100:150)]

BNf = as.data.frame(fNtt)
row.names(BNf) <- NULL 
BNf2 <- as.matrix(t(na.omit(t(BNf))))
# BNf$Phe <- as.numeric(BNf$Phe)
#knnf = MLearn( Phe ~ ., data=BNf, knnI(k=1,l=0),TrainInd)
knnf = MLearn( Phe ~ ., data=BNf2, randomForestI, TrainInd,ntree=600,na.action = na.fail)
knnf = MLearn( Phe ~ ., data=BNf2, knnI(k=3,l=2), TrainInd,ntree=600 )

confuMat(knnf)


#df<-data.frame(y=rnorm(10),x1=rnorm(10),x2=rnorm(10))
#lm(y~.,df)
library("fda.usc")
out1=classif.knn(Phe ~ ., data=BNf2,knn=3)
summary.classif(out1)
#    PREDICTION knn
mtest<-phoneme[["test"]][1:100]
gtest<-phoneme[["classtest"]][1:100]
pred1=predict.classif(out1,mtest)
table(pred1,gtest)




Traintt = rowttests(ALLfilt_bcrneg[, TrainInd], "mol.biol")
ordTT = order(abs(Traintt$statistic), decreasing=TRUE)
fNtt = featureNames(ALLfilt_bcrneg)[ordTT[1:50]]

BNf = ALLfilt_bcrneg[fNtt,]
knnf = MLearn( mol.biol ~ ., data=BNf, knnI(k=1,l=0),TrainInd)
confuMat(knnf)








