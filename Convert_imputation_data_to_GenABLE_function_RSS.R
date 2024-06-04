# Authoer: Jun Chen
# Date:    8/24/2013 
# Goal:    Convert the data format to 

setwd("/home/chenju/project_RSS/RSS/")

snp_60k <- read.table("/home/chenju/project_gs/data/SNP_60K_Sorted.txt",header=TRUE,sep="\t")
snp_60k_chrZ <- subset(snp_60k, Chr=="Z") # This is for sex checking.
snp_60k_chrW <- subset(snp_60k, Chr=="W")

# file1
temp <- read.table(file="/home/chenju/project_RSS/RSS/RSS_60k_beagle_haplotype.bgl.RSS_60k_beagle_haplotype.txt.phased",head=FALSE,sep=" ",check.names=F,colClasses="character")
# Line_SNP <- t(temp)

temp <- as.matrix(temp,nrow = dim(temp)[1], ncol=dim(temp)[2], byrow=TRUE,)
Mydata_temp <- convert.beagle.GenABLE(temp)

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


# Mydata_temp[10:16,1:15]
# temp[10:16,1:15]


# sex estimation
Sex_results <- sex.checking(snp_60k_chrZ,Line_SNP)
write.table(Sex_results,file="Sex_Est_Jun.txt",append=FALSE,quote=FALSE,sep=" ",row.names=TRUE,col.names=TRUE)

# 
# Function for sex checking.
sex.checking <- function(snp_60k_chrZ,Line_SNP)         # The inpute genotype data is from 
{

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
	return(Sex_results)
}
# end of sex estimation.



# Summary of SNP
N_row <- dim(Mydata_temp)[1]
N_col <- dim(Mydata_temp)[2]

# Mydata_temp_status <- substr(Mydata_temp[1,2:N_col],1,3)
# Phe <- Mydata_temp_status
# myData <- cbind(Phe,Mydata_temp)

sex <- read.table("/home/chenju/project_RSS/RSS/Sex_Est_Jun.txt",header=FALSE,sep=" ")
colnames(sex) <- c("id","sex")

Mydata_temp_status <- substr(Mydata_temp[1,2:N_col],1,3)
Phenotype <- t(rbind(Mydata_temp[1,2:N_col],as.numeric(as.factor(sex[,2]))-1,as.numeric(as.factor(Mydata_temp_status))-1))


colnames(Phenotype) <- c("id","sex","aff")
write.table(Phenotype,file="phe0.dat",append = FALSE, quote = FALSE,sep = "\t", row.names = FALSE, col.names = TRUE )

snp_60k <- read.table("/home/chenju/project_gs/data/SNP_60K_Sorted.txt",header=FALSE,sep="\t")   # need to change here.
colnames(snp_60k) <- c("name", "chr", "pos")
Geno_data <- cbind(snp_60k,Mydata_temp[,c(-1)])

setwd("/home/chenju/project_RSS/RSS/") # change here
write.table(Geno_data,file="gen0iwos_impute.dat",append = FALSE, quote = FALSE,sep = "\t", row.names = FALSE, col.names = FALSE )
# replace NA to 00

library(GenABEL)
convert.snp.illumina(inf="gen0iwos_impute.dat",out="gen0iwos_impute.raw",strand="+")

# start here.
setwd("/home/chenju/project_RSS/RSS/")
library(GenABEL)
df <- load.gwaa.data(phe="phe0.dat",gen="gen0iwos_impute.raw",force=TRUE)


# GWAS here.
descriptives.trait(df)
descriptives.marker(df)
# summary(df)
# head(summary(gtdata(df[(aff==1), ])))
 # P110
an0 <- qtscore(aff, df, trait="binomial")
results_Pvalue <- descriptives.scan(an0,top=57636)
write.table(results_Pvalue,file="gen0iwos_Pvalue.dat",append = FALSE, quote = FALSE,sep = "\t", row.names = TRUE, col.names = TRUE )



bad_SNP <- c("Gga_rs15930034","Gga_rs16234267","Gga_rs14393117","Gga_rs13772727","Gga_rs16683806")
sub_data1 <- subset(descriptives.scan(an0,top=57636), Chromosome != "Z")
num1 <- match(bad_SNP,rownames(sub_data1))
num2 <- which(num1>0)
sub_data <- sub_data1[-num1[num2],]
dim(sub_data)


pdf_name <- "Prediction_Model_for_RSS_PCA_no_chrZ_bad_SNP.pdf"
library(gtools)
pdf(pdf_name,width=8,height=8)
par(mfrow=c(2,3))


# PCA plot
# topnum <-504    # 57636,14,34,54,104,504
# topnum <- dim(sub_data)[1]  # topnum <- dim(sub_data)[1] if this is 60K.
topnum <-500  

num <- match(rownames(sub_data)[1:topnum],df@gtdata@snpnames)
# num <- match(rownames(descriptives.scan(an0,top=topnum)),df@gtdata@snpnames)
#num <- match(rownames(descriptives.scan(an0,top=topnum)[descriptives.scan(an0,top=topnum)$Chromosome !="Z",]),df@gtdata@snpnames)
length(num)
df2 <- df[,num[c(-1:-3,-11)]]   # delete the first three value

gkin <- ibs(df2, weight = "freq")
gkin[1:10,1:10]

cps.full <- cmdscale(as.dist(.5 - gkin), eig = T, k = 10)
names(cps.full)
cps <- cps.full$points
plot(cps[,1], cps[,2], pch = df2@phdata$aff,col=df2@phdata$aff+3,main="Top 500 SNPs")
# plot(cps[,1], cps[,2], pch = df2@phdata$aff,col=df2@phdata$aff+3,main="60K SNPs")
legend("topright",, legend=c("case","control"), pch = c(1,0),col=c(4,3))




dev.off()
# ***************************************************************************************************












# ***************************************************************************************************
# nj tree plot

setwd("/home/chenju/project_RSS/RSS") # change here

Mydata_temp_status <- substr(Mydata_temp[1,2:N_col],1,3)
Phe <- Mydata_temp_status

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


pdf_name <- "Prediction_Model_for_RSS_NJ_tree.pdf"
library(gtools)
pdf(pdf_name,width=20,height=20)

library(ape)
par(mfrow=c(2,2))

num <- match(rownames(descriptives.scan(an0,top=103)),as.character(Mydata_temp[,1]))
# MySNP <- Mydata_temp[,c(sort(num))]
MySNP <- Mydata_temp[num[-1:-3],-1]    # Delete the first three value.
rownames(MySNP) <- Mydata_temp[num[-1:-3],1]    # Delete the first three value.
colnames(MySNP) <- Mydata_temp[1,2:dim(Mydata_temp)[2]]    # Delete the first three value.

d <- dist(as.matrix(t(MySNP))) 
# tr <- nj(d)
tr <- njs(d)    # missing value
# ls(tr)  tr$tip.label  ?plot.phylo
X <- rep(c("red", "green"), each = 95)
plot(tr, "u",tip.col = X,show.node.label = TRUE,cex=0.5)
title("Top 100 SNPs")

dev.off()
















