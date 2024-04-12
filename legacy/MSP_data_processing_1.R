#enter directory location
directory <- "legacy/R Scripts New"
#file "protein-peptides.csv" exported from PEAKS label-free quantification
lfq_filename <- "legacy/protein-peptides-lfq.csv"
#file "protein-peptides.csv" exported from PEAKS identification
id_filename <- "legacy/protein-peptides-id.csv"
#specify mass spec samples, replicates should be denoted with the same numbers
sample <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)



##Optioanl, whether to fix random number generation
#set.seed(12345)
numsample <- length(sample)
library(data.table)
library(qvalue)
library(fitdistrplus)
library(truncnorm)
library(gtools)

#set working directory, all files should be in this directory

#setwd(directory)

#function "FDRcorrection", find peptides < 1%FDR
FDRcorrection <- function(lfq, id){
        df <- read.csv(lfq, sep = ",", header=TRUE, row.names = NULL, stringsAsFactors = FALSE)
        refdata <- read.csv(id, sep = ",", header=TRUE, row.names = NULL, stringsAsFactors = FALSE)
        rels <- c(refdata[,4])
        ls <- NULL
        for (i in 1:nrow(df)){
                ls <- c(ls, df[i,1:4], refdata[(which(rels == df[i,4])[1]),6:12], df[i,5:ncol(df)],which(rels == df[i,4])[1])
        }
        data <- as.matrix.data.frame(matrix(ls,ncol=(ncol(df)+8),byrow=TRUE))
        cnames1 <- colnames(df[1:4])
        cnames2 <- colnames(refdata[6:12])
        cnames3 <- colnames(df[5:ncol(df)])
        colnames(data) <- c(cnames1,cnames2,cnames3,"matchwhich")
        data<-data[complete.cases(data[ , ncol(data)]),]
}

#function "matfornorm", filter out Quality <= 0.3 and prepare the file for Normalyzer
matfornorm <- function (sample,df){
        data <- NULL
        ls <- NULL
        cnames <- NULL
        
        # Literally just filtering by quality score in the most confusing way. 
        for(i in 1:nrow(df)){
                if (df[i,14] > 0.3) {
                        ls <- c(ls,df[i,1:(17+numsample)])
                }
        }
        data <- as.matrix.data.frame(matrix(ls,ncol=(17+numsample),byrow=TRUE))
      
        cnames <- colnames(as.data.frame(df)[1:(17+numsample)])
    
        # Changing the column names, first 
        row0 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,sample)
        data <- rbind(row0,cnames,gsub("-", 0, matrix(gsub("-", 0, data[1:nrow(data),1:ncol(data)]),ncol = ncol(data))))
        colnames(data) <- cnames
        return (data)
}

data1 <- FDRcorrection(lfq_filename, id_filename)
data2 <- matfornorm(sample,data1)
data3 <- data2[2:nrow(data2),]

# Replacing 0 with NA
data3[data3 == 0] <- NA

#Adding a row with the sample groupings. This tells which ones are replicates. 
row0 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,sample)
data3 <- rbind(row0,data3)
write.csv(data1, file = "protein-peptides_FDRcorrected.csv", append = FALSE,  sep = ",", row.names = FALSE)
write.table(data2, file = "protein-peptides_for_Normalyzer.txt", append = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(data3, file = "protein-peptides_for_NormalyzerDE.txt", append = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
#go to https://urldefense.com/v3/__http://normalyzer.immunoprot.lth.se/normalize.php__;!!Mih3wA!GfyTNgvo0HrQ_DxH-T4HHknVIdHSrz7DO3BA2PJM2ujRE2PIsppu1hj7Wu22Q5-7V2xwohwfhx8VKaR5iQ$ , upload "protein-peptides_for_Normalyzer.txt" for normalization
#go to https://urldefense.com/v3/__http://quantitativeproteomics.org/normalyzerde__;!!Mih3wA!GfyTNgvo0HrQ_DxH-T4HHknVIdHSrz7DO3BA2PJM2ujRE2PIsppu1hj7Wu22Q5-7V2xwohwfhx_8aXVPCQ$ , upload "protein-peptides_for_NormalyzerDE.txt" for normalization
