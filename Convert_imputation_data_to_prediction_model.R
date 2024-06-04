# Authoer: Jun Chen
# Date:    8/24/2013 
# Goal:    Convert the data format to 

setwd("/home/chenju/project_RSS/RSS/")
#setwd("Z:/project_RSS/RSS/")

# snp_60k <- read.table("Z:/project_gs/data/SNP_60K_Sorted.txt",header=TRUE,sep="\t")
snp_60k <- read.table("/home/chenju/project_gs/data/SNP_60K_Sorted.txt",header=TRUE,sep="\t")
snp_60k_chrZ <- subset(snp_60k, Chr=="Z") # This is for sex checking.
snp_60k_chrW <- subset(snp_60k, Chr=="W")

# file1
temp <- read.table(file="/home/chenju/project_RSS/RSS/RSS_60k_beagle_haplotype.bgl.RSS_60k_beagle_haplotype.txt.phased",head=FALSE,sep=" ",check.names=F,colClasses="character")
# temp <- read.table(file="Z:/project_RSS/RSS/RSS_60k_beagle_haplotype.bgl.RSS_60k_beagle_haplotype.txt.phased",head=FALSE,sep=" ",check.names=F,colClasses="character")
# Line_SNP <- t(temp)

convert.beagle.GenABLE <- function(d)         # The inpute genotype data is from 
{
	N_row <- dim(d)[1]
	N_col <- dim(d)[2]
	hap.ori <- array(0,c(N_row,N_col/2))
	hap.ori[,]=paste(d[,seq(1,N_col,by=2)],d[,seq(2,N_col,by=2)],sep="")
	hap.ori[1,] <- d[1,seq(2,N_col,by=2)]
	hap.ori[,1] <- d[,2]
	return(hap.ori)
}

temp <- as.matrix(temp,nrow = dim(temp)[1], ncol=dim(temp)[2], byrow=TRUE,)
Mydata_temp <- convert.beagle.GenABLE(temp)


# Summary of SNP

# GWAS results from GenABLE. (The method is based on Plink)
results_Pvalue <- read.table(file="gen0iwos_Pvalue.dat",head=TRUE,sep = "\t", check.names=F,colClasses="character",row.names=1)
# pickup_top 50

#source("http://bioconductor.org/biocLite.R")
#biocLite("genefilter")

library("Biobase")
library("genefilter")
library("bioDist")
library("MLInterfaces")
library("randomForest")
library("MASS")

genotype11<-"AA"
genotype12<-"AB"
genotype21<-"BA"
genotype22<-"BB"
for(i in 1:dim(Mydata_temp)[2])
{
  if(length(which(Mydata_temp[,i]==as.character(genotype11)))>0) Mydata_temp[c(which(Mydata_temp[,i]==as.character(genotype11))),i]="0"
  if(length(which(Mydata_temp[,i]==as.character(genotype12)))>0) Mydata_temp[c(which(Mydata_temp[,i]==as.character(genotype12))),i]="1"
  if(length(which(Mydata_temp[,i]==as.character(genotype21)))>0) Mydata_temp[c(which(Mydata_temp[,i]==as.character(genotype21))),i]="1"
  if(length(which(Mydata_temp[,i]==as.character(genotype22)))>0) Mydata_temp[c(which(Mydata_temp[,i]==as.character(genotype22))),i]="2"
}

# rownames(MySNP) <- rownames(myData)

for(j in 2:dim(Mydata_temp)[2])
{
  Mydata_temp[2:dim(Mydata_temp)[1],j] <- as.numeric(Mydata_temp[2:dim(Mydata_temp)[1],j])
}

N_row <- dim(Mydata_temp)[1]
N_col <- dim(Mydata_temp)[2]
Mydata_temp_status <- substr(Mydata_temp[1,2:N_col],1,3)
Phe <- Mydata_temp_status

Mydata <- t(Mydata_temp)[-1,-1]
rownames(Mydata) <- Mydata_temp[1,2:N_col]
colnames(Mydata) <- Mydata_temp[2:N_row,1]

# cluster here, see what happened.
# Traintt = rowttests(Mydata[,TrainInd],f2[TrainInd])
# a <- table(Phe[TrainInd],Mydata[i,TrainInd])
# ordTT = order(abs(Traintt$statistic), decreasing=TRUE)

Mydata2 <- cbind(Phe,Mydata)
num <- match(rownames(results_Pvalue)[4:53],colnames(Mydata2))
RSS = as.data.frame(Mydata2[,c(1,num)])


# we can do a loop here, for 100 simulation.
num_train <- 50
neg = which(Phe == "neg")
pos = which(Phe == "pos")
S1 = sample(neg, num_train, replace=FALSE)
S2 = sample(pos, num_train, replace = FALSE)
TrainInd = c(S1, S2)
TestInd = setdiff(1:190, TrainInd)
# #f2 <- factor(Phe)


knnf = MLearn( Phe ~ ., data=RSS, knnI(k=1,l=0),TrainInd)
confuMat(knnf)
knnf = MLearn( Phe ~ ., data=RSS, randomForestI, TrainInd,ntree=600)
confuMat(knnf) 
knnf <- MLearn(Phe ~ ., data=RSS, svmI, TrainInd)
confuMat(knnf) 
knnf = MLearn( Phe ~ ., data=RSS, naiveBayesI,TrainInd)
confuMat(knnf)  # on held-out test data
knnf = MLearn( Phe ~ ., data=RSS, rpartI, TrainInd)
confuMat(knnf)  # on held-out test data




# Method 1
knnf = MLearn( Phe ~ ., data=RSS, knnI(k=1,l=0),TrainInd)
confuMat(knnf)
# RObject(knnf)


pdf_name <- "Prediction_Model_for_RSS_randomForest.pdf"
library(gtools)
pdf(pdf_name,width=20,height=20)

# Method 2. randomForest methods.   It is great.
knnf = MLearn( Phe ~ ., data=RSS, randomForestI, TrainInd,ntree=600)
confuMat(knnf) 
RObject(knnf)
names(RObject(knnf))

dev.off()

## examples for predict
library("e1071")
knnf <- MLearn(Phe ~ ., data=RSS, svmI, TrainInd)
confuMat(knnf) 
predict(knnf, RSS[TestInd,])
results <- predict(knnf, RSS[TestInd,-1])
results$testPredictions
results$testScores





pdf_name <- "Prediction_Model_for_RSS_support_vector_machine.pdf"
library(gtools)
pdf(pdf_name,width=10,height=10)
matplot(results$testScores,)
dev.off()


# add plot forr clust pca and pca



# ## examples for predict
# clout <- MLearn(type~., sample.ExpressionSet[100:250,], svmI , 1:16)
# predict(clout, sample.ExpressionSet[100:250,17:26])


# ------------
# Method 3. 
knnf = MLearn( Phe ~ ., data=RSS, dldaI, TrainInd)
confuMat(knnf) 

knnf = MLearn( Phe ~ ., data=RSS, ldaI, TrainInd)
confuMat(knnf) 

knnf = MLearn( Phe ~ ., data=RSS, ldaI.predParms(method="debiased"), TrainInd)
confuMat(knnf) 
 
library("ipred")
knnf = MLearn( Phe ~ ., data=RSS, sldaI, TrainInd)    # need package ipred
confuMat(knnf) 

# knnf = MLearn( Phe ~ ., data=RSS, qdaI, TrainInd)    # need package ipred
# confuMat(knnf) 

# Method 6   support vector machine
library("e1071")
knnf = MLearn( Phe ~ ., data=RSS, svmI, TrainInd)    # need package e1071
confuMat(knnf)  # on held-out test data

# Method 4
knnf = MLearn( Phe ~ ., data=RSS, rpartI, TrainInd)
# confuMatTrain(knnf)  # on training data
confuMat(knnf)  # on held-out test data

# Method 5
knnf = MLearn(Phe ~ ., data=RSS, glmI.logistic(threshold=0.5),TrainInd,family=binomial)
confuMat(knnf)  # on held-out test data


# Method 7
## recode data for RAB
#nsp = ifelse(crabs$sp=="O", -1, 1)
#nsp = factor(nsp)
#ncrabs = cbind(nsp,crabs)
#rab1 = MLearn(nsp~CW+RW, data=ncrabs, RABI, kp, maxiter=10)
#rab1
#
# knnf = MLearn( Phe ~ ., data=RSS, RABI, TrainInd,maxiter=20, maxdepth=2)  # need recode the data
# confuMat(knnf)  # on held-out test data

# Method 8
library(ada)
knnf = MLearn( Phe ~ ., data=RSS, adaI, TrainInd,)   # need packages ada
confuMat(knnf)  # on held-out test data


# pdf_name <- "Prediction_Model_for_RSS_randomForest.pdf"
# library(gtools)
# pdf(pdf_name,width=8,height=8)

# library(ada)
# knnf = MLearn( Phe ~ ., data=RSS, rdacvI, TrainInd)   # need packages ada
# confuMat(knnf)  # on held-out test data
# RObject(knnf)
# plotXvalRDA(knnf)  # special interface to plots of parameter space

# dev.off()


# Method 8
knnf = MLearn( Phe ~ ., data=RSS, lvqI, TrainInd)
confuMat(knnf)  # on held-out test data

# Method 8
knnf = MLearn( Phe ~ ., data=RSS, naiveBayesI,TrainInd)
confuMat(knnf)  # on held-out test data

# Method 8
knnf = MLearn( Phe ~ ., data=RSS, baggingI, TrainInd)
confuMat(knnf)  # on held-out test data


# new mboost interface -- you MUST supply family for nonGaussian response
#
# require(party)  # trafo ... killing cmd check
# blb.1 = MLearn(sp~CW+RW+FL, data=crabs, blackboostI, kp, family=mboost::Binomial() )
# confuMat(blb.1)
















