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
# results_Pvalue <- read.table(file="gen0iwos_Pvalue.dat",head=TRUE,sep = "\t", check.names=F,colClasses="character",row.names=1)
# # pickup_top 50

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


# bad_marker <- c("Gga_rs15930034","Gga_rs16234267","Gga_rs14393117","Gga_rs13772727","Gga_rs16683806")
# num_bad <- match(bad_marker,rownames(results_Pvalue))

setwd("/home/chenju/project_RSS/RSS/")
library(GenABEL)
df2 <- load.gwaa.data(phe="phe0.dat",gen="gen0iwos_impute.raw",force=TRUE)


# simulation start here.
Accuracy <- array("NA",c(50*5*3,3))
c<-0
for (i in 1:50)
{

# we can do a loop here, for 100 simulation.
num_train <- 70
neg = which(Phe == "neg")
pos = which(Phe == "pos")
S1 = sample(neg, num_train, replace=FALSE)
S2 = sample(pos, num_train, replace = FALSE)
TrainInd = c(S1, S2)
TestInd = setdiff(1:190, TrainInd)
# #f2 <- factor(Phe)

# This is GWAS
df <- df2[TrainInd,]

an0 <- qtscore(aff, df, trait="binomial")

bad_SNP <- c("Gga_rs15930034","Gga_rs16234267","Gga_rs14393117","Gga_rs13772727","Gga_rs16683806")
sub_data1 <- subset(descriptives.scan(an0,top=57636), Chromosome != "Z")
num1 <- match(bad_SNP,rownames(sub_data1))
num2 <- which(num1>0)
sub_data <- sub_data1[-num1[num2],]
dim(sub_data)


for (j in 1:3)
{
# top 50
if(j==1) Topnum <-10
# if(j==1) Topnum <-500   # there are some bad SNP in there, It is not good to predict.
if(j==2) Topnum <-50
if(j==3) Topnum <-100
# Topnum <- 50
num <- match(rownames(sub_data)[1:Topnum],colnames(Mydata2))
RSS = as.data.frame(Mydata2[,c(1,num)])


c<-try(MLearn( Phe ~ ., data=RSS, randomForestI, TrainInd,ntree=600)@testScores,T)
if(any(nchar(c)>40)){next}

# prediction model.  start here
knnf = MLearn( Phe ~ ., data=RSS, randomForestI, TrainInd,ntree=600)  # random Forest
# confuMat(knnf)
Accuracy[15*(i-1)+j*5-4,3]=paste("Top ",Topnum," SNPs",sep="")
Accuracy[15*(i-1)+j*5-4,1]="Random Forest"
Accuracy[15*(i-1)+j*5-4,2]=(confuMat(knnf)[1,2]+confuMat(knnf)[2,1])/(confuMat(knnf)[1,1]+confuMat(knnf)[1,2]+confuMat(knnf)[2,1]+confuMat(knnf)[2,2])
knnf <- MLearn(Phe ~ ., data=RSS, svmI, TrainInd)    # SVM
# confuMat(knnf) 
Accuracy[15*(i-1)+j*5-3,3]=paste("Top ",Topnum," SNPs",sep="")
Accuracy[15*(i-1)+j*5-3,1]="SVM"
Accuracy[15*(i-1)+j*5-3,2]=(confuMat(knnf)[1,2]+confuMat(knnf)[2,1])/(confuMat(knnf)[1,1]+confuMat(knnf)[1,2]+confuMat(knnf)[2,1]+confuMat(knnf)[2,2])
knnf = MLearn( Phe ~ ., data=RSS, naiveBayesI,TrainInd)  # naiveBayes
# confuMat(knnf)  # on held-out test data
Accuracy[15*(i-1)+j*5-2,3]=paste("Top ",Topnum," SNPs",sep="")
Accuracy[15*(i-1)+j*5-2,1]="Naive Bayes"
Accuracy[15*(i-1)+j*5-2,2]=(confuMat(knnf)[1,2]+confuMat(knnf)[2,1])/(confuMat(knnf)[1,1]+confuMat(knnf)[1,2]+confuMat(knnf)[2,1]+confuMat(knnf)[2,2])
knnf = MLearn( Phe ~ ., data=RSS, knnI(k=1,l=0),TrainInd)  # knnI
# confuMat(knnf)
Accuracy[15*(i-1)+j*5-1,3]=paste("Top ",Topnum," SNPs",sep="")
Accuracy[15*(i-1)+j*5-1,1]="KNN"
Accuracy[15*(i-1)+j*5-1,2]=(confuMat(knnf)[1,2]+confuMat(knnf)[2,1])/(confuMat(knnf)[1,1]+confuMat(knnf)[1,2]+confuMat(knnf)[2,1]+confuMat(knnf)[2,2])
knnf = MLearn( Phe ~ ., data=RSS, rpartI, TrainInd) # rpart
# confuMat(knnf)  # on held-out test data
Accuracy[15*(i-1)+j*5,3]=paste("Top ",Topnum," SNPs",sep="")
Accuracy[15*(i-1)+j*5,1]="RPRT"
Accuracy[15*(i-1)+j*5,2]=(confuMat(knnf)[1,2]+confuMat(knnf)[2,1])/(confuMat(knnf)[1,1]+confuMat(knnf)[1,2]+confuMat(knnf)[2,1]+confuMat(knnf)[2,2])

}
print(i)
}



# some SNP only have two genotype, that is not good.

Accuracy=as.data.frame(Accuracy)
colnames(Accuracy)=c("Method","Acc","SNP")
# Accuracy2=subset(Accuracy, !is.na(Method) | !is.na(Acc))

Accuracy2=Accuracy[which(Accuracy$Method != "NA"),]
dim(Accuracy2)/15

Accuracy2$Method=factor(Accuracy2$Method)
Accuracy2$Acc=as.numeric(unlist(as.character(Accuracy2$Acc)))*100
Accuracy2$SNP=factor(Accuracy2$SNP)


levels(Accuracy2$SNP) <- c("Top 10 SNPs","Top 50 SNPs","Top 100 SNPs")
# levels(Accuracy2$SNP) <- c("Top 50 SNPs","Top 100 SNPs","Top 500 SNPs")
# plot 1
library(ggplot2)
library(lattice)

pdf_name <- "Boxplot_for_RSS_5_prediction_methods_70_top_multi_SNPs_lattice.pdf"
pdf(pdf_name,width=8,height=8)
# boxplot(Acc~Method,data=Accuracy2, main="Prediction by 50 simulations",col=(c("#D7191C","#FDAE61","#FFFFBF","#ABDDA4","#2B83BA")),na.action = NULL,
   # xlab="Method", ylab="Error rate (%)")  
bwplot(Acc ~ Method | factor(SNP),Accuracy2, varwidth = TRUE, layout = c(1, 3), groups=Method,xlab="Methods", ylab="Error rate (%)")
bwplot(Acc ~ Method | factor(SNP),Accuracy2, varwidth = TRUE, layout = c(3, 1), groups=Method,xlab="Methods", ylab="Error rate (%)")
dev.off()   



pdf_name <- "Boxplot_for_RSS_5_prediction_methods_70_top_multi_SNPs_ggplot2.pdf"
pdf(pdf_name,width=12,height=8)
 ggplot(Accuracy2, aes(factor(Method), Acc,fill = factor(SNP))) + geom_boxplot() +labs(title = "Prediction by 50 simulations",x="Methods",y="Error rate (%)")  
ggplot(Accuracy2, aes(Method, Acc,fill = SNP)) + geom_boxplot() + facet_wrap(~ SNP, ncol = 1)+labs(title = "Prediction by 50 simulations",x="Methods",y="Error rate (%)")   
ggplot(Accuracy2, aes(factor(Method), Acc,fill = SNP)) + geom_boxplot() + facet_grid(~SNP) + labs(title = "Prediction by 50 simulations",x="Methods",y="Error rate (%)") 
dev.off()   


# lattice

# > pl <- bwplot(Days ~ log(FSC.H), data = gvhd10, xlab = "log(Forward Scatter)",
# +     ylab = "Days Past Transplant")
# > print(pl)
# ggplot2
# > p <- ggplot(gvhd10, aes(factor(Days), log(FSC.H)))
# > pg <- p + geom_boxplot() + coord_flip() + labs(y = "log(Forward Scatter)",
# +     x = "Days Past Transplant")
# > print(pg)
# lattice
# > pl <- bwplot(gcsescore^2.34 ~ gender | factor(score),
# +     Chem97, varwidth = TRUE, layout = c(6, 1), ylab = "Transformed GCSE score")
# > print(pl)
# ggplot2
# > p <- ggplot(Chem97, aes(factor(gender), gcsescore^2.34))
# > pg <- p + geom_boxplot() + facet_grid(~score) + ylab("Transformed GCSE score")
# > print(pg)
# lattice
# > pl <- bwplot(factor(score) ~ gcsescore | gender, data = Chem97,
# +     xlab = "Average GCSE Score")
# > print(pl)
# ggplot2
# > pg <- ggplot(Chem97, aes(factor(score), gcsescore)) +
# +     geom_boxplot() + coord_flip() + ylab("Average GCSE score") +
# +     facet_wrap(~gender)
# > print(pg)
# # http://learnr.wordpress.com/2009/06/29/ggplot2-version-of-figures-in-lattice-multivariate-data-visualization-with-r-part-2/
# > pg <- ggplot(barley, aes(yield, variety, colour = year)) +
# +     geom_point() + facet_wrap(~site, ncol = 1)
# > print(pg)


# library(RColorBrewer)
# a=brewer.pal(5,"Spectral")
# > a
# [1] "#D7191C" "#FDAE61" "#FFFFBF" "#ABDDA4" "#2B83BA"




  
  
  
  
pdf_name <- "Prediction_Model_for_RSS_support_vector_machine_70.pdf"
library(gtools)
pdf(pdf_name,width=7,height=7)
# matplot(results$testScores,pch=21:22,col=2:3,bg=2:3)
matplot(knnf$testScores,col=2:3,xlab="Individuals",ylab="Predict disease status")
dev.off()









# # This is for plot

# # machine learning methods. support vector machine.
# ## examples for predict   # running here.
# library("e1071")
# knnf <- MLearn(Phe ~ ., data=RSS, svmI, TrainInd)
# confuMat(knnf) 
# # predict(knnf, RSS[TestInd,])
# results <- predict(knnf, RSS[TestInd,-1])
# results$testPredictions
# results$testScores

# pdf_name <- "Prediction_Model_for_RSS_support_vector_machine_70.pdf"
# library(gtools)
# pdf(pdf_name,width=7,height=7)
# # matplot(results$testScores,pch=21:22,col=2:3,bg=2:3)
# matplot(results$testScores,col=2:3,xlab="Individuals",ylab="Predict disease status")
# dev.off()

# # 


# # Method 1
# knnf = MLearn( Phe ~ ., data=RSS, knnI(k=1,l=0),TrainInd)
# confuMat(knnf)
# # RObject(knnf)


# pdf_name <- "Prediction_Model_for_RSS_randomForest.pdf"
# library(gtools)
# pdf(pdf_name,width=20,height=20)

# # Method 2. randomForest methods.   It is great.
# knnf = MLearn( Phe ~ ., data=RSS, randomForestI, TrainInd,ntree=600)
# confuMat(knnf) 
# RObject(knnf)
# names(RObject(knnf))

# dev.off()



# # add plot forr clust pca and pca



# # ## examples for predict
# # clout <- MLearn(type~., sample.ExpressionSet[100:250,], svmI , 1:16)
# # predict(clout, sample.ExpressionSet[100:250,17:26])


# # ------------
# # Method 3. 
# knnf = MLearn( Phe ~ ., data=RSS, dldaI, TrainInd)
# confuMat(knnf) 

# knnf = MLearn( Phe ~ ., data=RSS, ldaI, TrainInd)
# confuMat(knnf) 

# knnf = MLearn( Phe ~ ., data=RSS, ldaI.predParms(method="debiased"), TrainInd)
# confuMat(knnf) 
 
# library("ipred")
# knnf = MLearn( Phe ~ ., data=RSS, sldaI, TrainInd)    # need package ipred
# confuMat(knnf) 

# # knnf = MLearn( Phe ~ ., data=RSS, qdaI, TrainInd)    # need package ipred
# # confuMat(knnf) 

# # Method 6   support vector machine
# library("e1071")
# knnf = MLearn( Phe ~ ., data=RSS, svmI, TrainInd)    # need package e1071
# confuMat(knnf)  # on held-out test data

# # Method 4
# knnf = MLearn( Phe ~ ., data=RSS, rpartI, TrainInd)
# # confuMatTrain(knnf)  # on training data
# confuMat(knnf)  # on held-out test data

# # Method 5
# knnf = MLearn(Phe ~ ., data=RSS, glmI.logistic(threshold=0.5),TrainInd,family=binomial)
# confuMat(knnf)  # on held-out test data


# # Method 7
# ## recode data for RAB
# #nsp = ifelse(crabs$sp=="O", -1, 1)
# #nsp = factor(nsp)
# #ncrabs = cbind(nsp,crabs)
# #rab1 = MLearn(nsp~CW+RW, data=ncrabs, RABI, kp, maxiter=10)
# #rab1
# #
# # knnf = MLearn( Phe ~ ., data=RSS, RABI, TrainInd,maxiter=20, maxdepth=2)  # need recode the data
# # confuMat(knnf)  # on held-out test data

# # Method 8
# library(ada)
# knnf = MLearn( Phe ~ ., data=RSS, adaI, TrainInd,)   # need packages ada
# confuMat(knnf)  # on held-out test data


# # pdf_name <- "Prediction_Model_for_RSS_randomForest.pdf"
# # library(gtools)
# # pdf(pdf_name,width=8,height=8)

# # library(ada)
# # knnf = MLearn( Phe ~ ., data=RSS, rdacvI, TrainInd)   # need packages ada
# # confuMat(knnf)  # on held-out test data
# # RObject(knnf)
# # plotXvalRDA(knnf)  # special interface to plots of parameter space

# # dev.off()


# # Method 8
# knnf = MLearn( Phe ~ ., data=RSS, lvqI, TrainInd)
# confuMat(knnf)  # on held-out test data

# # Method 8
# knnf = MLearn( Phe ~ ., data=RSS, naiveBayesI,TrainInd)
# confuMat(knnf)  # on held-out test data

# # Method 8
# knnf = MLearn( Phe ~ ., data=RSS, baggingI, TrainInd)
# confuMat(knnf)  # on held-out test data


# # new mboost interface -- you MUST supply family for nonGaussian response
# #
# # require(party)  # trafo ... killing cmd check
# # blb.1 = MLearn(sp~CW+RW+FL, data=crabs, blackboostI, kp, family=mboost::Binomial() )
# # confuMat(blb.1)
















