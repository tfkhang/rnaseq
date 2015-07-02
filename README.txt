
R code documentation
Source: Khang TF, Lau CY. (2015) Getting the most out of RNA-seq data analysis. PeerJ PrePrints 3:e1463 https://dx.doi.org/10.7287/peerj.preprints.1198v1
Author: Ching Yee Lau
Version: 1.0
Date: 26 June 2016 


The 'DEA' folder contains the R implementation of Differential Expression Analysis (DEA) of raw RNA-Seq count data using DESeq, DESeq2, NOISeq, edgeR, ASC and simple Z-test. For a comprehensive description of this function, please refer to the DEA_manual.txt in 'DEA' folder.  

'GFOLD.R' file contains the pipeline for implementing GFOLD analysis in Khang & Lau (2015)

Please note that the DEA function and 'GFOLD.R' are for running a single sample combination of simulated data in Khang & Lau (2015).

For running all the simulated data, one might consider using 'apply', 'lapply', 'sapply', 'mapply' function or 'for' loop in R.    

Besides the simulations done on the Cheung and Bottomly data sets in Khang & Lau (2015), DEA function and 'GFOLD.R' can be applied to any other raw RNA-Seq count data sets.


