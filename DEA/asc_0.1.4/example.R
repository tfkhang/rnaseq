source("./asc_0.1.4.R")
load("exampleData.rda")
### The example data set has two long vectors, X1 and X2, each of 20,000 genes
### X1 and X2 must of the same length and should represent genes in the same order
### To estimate the posterior mean of log fold change 
D <- DGE(X1,X2)
### This function returns the estimated log fold change, as well as the hyperparamters
###
### identify a subset of genes of interest.
### For example, those with log fold change greater than 0.5
k1=which(abs(D$delta[,1])>.5)
D$delta[k1] 
### obtain estimated probability that the fold chnage is greater than
### a cut off. Here the example is 2 (notice this is fold change, not fold change).
postp1=PostProb(D$delta[k1,],X1[k1],X2[k1],D$pars,d0=2)
round(postp1,3)
###ANother example
### find the top 10 genes with the greatest estimated fold change
k2=order(-abs(D$delta[,1]))[1:10]
### estimate posterior probability that the fold change is more than 50%
postp2=PostProb(D$delta[k2,],X1[k2],X2[k2],D$pars,d0=1.5)
round(postp2,3)


