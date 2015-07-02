################################################################################################################
# R code documentation                                                                                         #
# Khang TF, Lau CY. (2015) Getting the most out of RNA-seq 
# data analysis. PeerJ PrePrints 3:e1463 
# https://dx.doi.org/10.7287/peerj.preprints.1198v1
# Author: Ching Yee Lau                                                                                        #
# Version: 1.0                                                                                                 #
# Date: 26 June 2016                                                                                           #
################################################################################################################


## Pipeline for implementing GFOLD in Khang & Lau (2015)

# Load DESeq package

library(DESeq)

# Raw RNA-Seq count data
## Source example: cheung.eset (ExpressionSet object) from the ReCount website 'http://bowtie-bio.sourceforge.net/recount/' 
## For demo, only use subset of cheung.eset with 3 female samples and 3 male samples
## To analyse other data sets, simply replace the "countData" and "conditions" objects below.
## For the details of "countData" and "conditions", see the newCountDataSet documentation ('?newCountDataSet').

countData <- exprs(cheung.eset)[,c(1:3,18:20)]
conditions <- cheung.eset$gender[c(1:3,18:20)]


# Perform normalization using DESeq method

cds <- newCountDataSet(countData,conditions)
cds <- estimateSizeFactors(cds)


# Omit the genes with zero counts across all the samples

cds@featureData@data <- cds@featureData@data[names(which(rowSums(counts(cds,normalized=TRUE))!=0)),]
counts(cds) <- counts(cds)[rownames(cds@featureData@data),]


# Create text files for use as input in GFOLD analysis
## Each sample has one text file named as "sample<num>.txt", where <num> is the column number in 'countData'
## For the example here, six text files were created: 'sample1.txt','sample2.txt', ..., 'sample6.txt', representing column 1 to 6 of 'countData' respectively

mapply(function(x){write.table(cbind(NA,counts(cds,normalized=TRUE)[,x],NA,NA),file=paste("sample",x,".txt",sep=""),sep="\t",col.names=FALSE)},x=1:ncol(cds))


## The open source C/C++ program for GFOLD analysis is available at 'http://www.tongji.edu.cn/∼zhanglab/GFOLD/index.html'
## After correctly installing GFOLD, run the following command at the Linux terminal for GFOLD analysis of the example data here: 

# > gfold diff -s1 sample1,sample2,sample3 -s2 sample4,sample5,sample6 -suf .txt -o 123vs456.diff -norm NO

## Note that the above command '-norm NO' stands for no normalization since the input files are normalized data for each sample
## Using other data sets, the above command may need to be modified and the output filename can be changed to any preferred name
## For a comprehensive description of GFOLD analysis, please refer to the GFOLD manual available at the website mentioned above


# Read output file into R

GFOLD.out <- read.delim("123vs456.diff",row.names="X.GeneSymbol",skip=ncol(cds)+7)

## If different filename is used, please change "123vs456.diff" to the corresponse filename

# Get the differentially expressed genes (DEG)
## The selection criterion is |GFOLD value| >= 1 , this cutoff can be changed 

GFOLD.deg <- rownames(GFOLD.out)[abs(GFOLD.out$GFOLD.0.01.)>=1]


