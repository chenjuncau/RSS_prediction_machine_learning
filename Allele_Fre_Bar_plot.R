
# Jun Chen
# Aug 9, 2013


setwd("C:/JChen/Project_RSS/Results")
SNP_M_F <- read.table(file="RSS_M_F_Chi2_top100.txt",head=FALSE,sep=" ")
colnames(SNP_M_F) <- c("chrom","pos","ref","alt","QV","RD","QS","Name","altMAF","Afre","ref1","alt1","SNP1","SNP2","ref2","alt2","oddsratio","P_fisher","P_chi2","chrom_M","pos_M","ref_M","alt_M","QV_M","RD_M","QS_M","Name_M","altMAF_M","Afre_M","ref1_M","alt1_M","SNP1_M","SNP2_M","ref2_M","alt2_M","oddsratio_M","P_fisher_M","P_chi2_M")

# SNP_M_F[1:10,]

# Result <- array(0,c(2,10))
# for(j in 1:10)
# {
# Result[1,j] <- SNP_M_F$alt1/(SNP_M_F$alt1+SNP_M_F$ref1)
# Result[2,j] <- SNP_M_F$alt2/(SNP_M_F$alt2+SNP_M_F$ref2)
# }


pdf("SNP_Fre_20_plot.pdf",width=20,height=8)
par(mfrow=c(2,1))
# Female pool

Allele_1 <- SNP_M_F$alt1/(SNP_M_F$alt1+SNP_M_F$ref1)
Allele_2 <- SNP_M_F$ref1/(SNP_M_F$alt1+SNP_M_F$ref1)

Result <- array(0,c(2,20))
Result[1,] <- Allele_1[1:20]
Result[2,] <- Allele_2[1:20]
colnames(Result) <- paste("SNP",1:20,sep="")

colors <- c("black", "red")
barplot(height = Result,beside = FALSE,col = colors,main="Positive Allele Frequency for Female Pool",cex.lab=0.8)#,legend.text = c("A", "B"))#,args.legend ="topright")
# legend("topright", inset=c(-1,0), legend = c("A", "B"), fill = colors, bty = "n")


Allele_1 <- SNP_M_F$alt2/(SNP_M_F$alt2+SNP_M_F$ref2)
Allele_2 <- SNP_M_F$ref2/(SNP_M_F$alt2+SNP_M_F$ref2)

Result <- array(0,c(2,20))
Result[1,] <- Allele_1[1:20]
Result[2,] <- Allele_2[1:20]
colnames(Result) <- paste("SNP",1:20,sep="")

colors <- c("black", "red")
barplot(height = Result,beside = FALSE,col = colors,main="Negative Allele Frequency for Female Pool",cex.lab=0.8)#,legend.text = c("A", "B"))#,args.legend ="topright")
# legend("topright", inset=c(-1,0), legend = c("A", "B"), fill = colors, bty = "n")

		
# Male		


Allele_1 <- SNP_M_F$alt1_M/(SNP_M_F$alt1_M+SNP_M_F$ref1_M)
Allele_2 <- SNP_M_F$ref1_M/(SNP_M_F$alt1_M+SNP_M_F$ref1_M)

Result <- array(0,c(2,20))
Result[1,] <- Allele_1[1:20]
Result[2,] <- Allele_2[1:20]
colnames(Result) <- paste("SNP",1:20,sep="")

colors <- c("black", "red")
barplot(height = Result,beside = FALSE,col = colors,main="Positive Allele Frequency for Male Pool",cex.lab=0.8)#,legend.text = c("A", "B"))#,args.legend ="topright")
# legend("topright", inset=c(-1,0), legend = c("A", "B"), fill = colors, bty = "n")


Allele_1 <- SNP_M_F$alt2_M/(SNP_M_F$alt2_M+SNP_M_F$ref2_M)
Allele_2 <- SNP_M_F$ref2_M/(SNP_M_F$alt2_M+SNP_M_F$ref2_M)

Result <- array(0,c(2,20))
Result[1,] <- Allele_1[1:20]
Result[2,] <- Allele_2[1:20]
colnames(Result) <- paste("SNP",1:20,sep="")

colors <- c("black", "red")
barplot(height = Result,beside = FALSE,col = colors,main="Negative Allele Frequency for Male Pool",cex.lab=0.8)#,legend.text = c("A", "B"))#,args.legend ="topright")
# legend("topright", inset=c(-1,0), legend = c("A", "B"), fill = colors, bty = "n")

dev.off()
