
library(DESeq)
library(DESeq2)
library(NOISeq)
library(edgeR)
source("DEA/asc_0.1.4/asc_0.1.4.R")
source("DEA/ztest.R")


DEA <- function(countData,conditions,replicates=NULL,bcv){

if(is.null(replicates)) stop ('need to specify if experiment has replicates (TRUE) or not (FALSE)')

# perform normalization by DESeq

cds <- newCountDataSet(countData,conditions)
cds <- estimateSizeFactors(cds)


# omit the genes which have no expression (value=0) across all the samples

cds@featureData@data <- cds@featureData@data[names(which(rowSums(counts(cds,normalized=TRUE))!=0)),]
counts(cds) <- counts(cds)[rownames(cds@featureData@data),]


# perform DESeq to indentify DE genes

if(replicates){

cds <- estimateDispersions(cds)

}

else{

cds <- estimateDispersions(cds,method="blind",sharingMode="fit-only")

}

deseq.out <- nbinomTest(cds,levels(conditions(cds))[1],levels(conditions(cds))[2])

deseq.deg <- deseq.out[-log10(deseq.out$pval)>=2/abs(deseq.out$log2FoldChange),"id"]


# perform DESeq2 to indentify DE genes

deseq2.out <- results(DESeq(DESeqDataSetFromMatrix(countData=counts(cds),colData=pData(cds),design=~condition)))

deseq2.deg <- rownames(deseq2.out)[-log10(deseq2.out$pvalue)>=2/abs(deseq2.out$log2FoldChange)]


# perform NOISeq to indentify DE genes

if(replicates){

noiseq.out <- noiseqbio(readData(counts(cds,normalized=TRUE),as.data.frame(conditions(cds))),norm="n",factor="conditions(cds)",filter=0)

noiseq.deg <- rownames(degenes(noiseq.out,q=0.95))

}

else{

noiseq.out <- noiseq(readData(counts(cds,normalized=TRUE),as.data.frame(conditions(cds))),norm="n",replicates="no",factor="conditions(cds)")

noiseq.deg <- rownames(degenes(noiseq.out,q=0.9))

}


# perform edgeR to indentify DE genes

if(replicates){

edger.out <- topTags(exactTest(estimateTagwiseDisp(estimateCommonDisp(DGEList(counts=counts(cds,normalized=TRUE),group=conditions(cds))))),n=Inf,adjust.method="none")

}

else{

edger.out <- topTags(exactTest(DGEList(counts=counts(cds,normalized=TRUE),group=conditions(cds)),dispersion=bcv^2),n=Inf,adjust.method="none")

}

edger.deg <- rownames(edger.out$table)[-log10(edger.out$table$PValue)>=2/abs(edger.out$table$logFC)]


# perform ASC to indentify DE genes

if(replicates){

asc.deg <- NA

}

else{

asc.out <- try(DGE(round(counts(cds,normalized=TRUE)[,1]),round(counts(cds,normalized=TRUE)[,2])),T)

asc.out$postp_d0_2 <- try(PostProb(asc.out$delta,round(counts(cds,normalized=TRUE)[,1]),round(counts(cds,normalized=TRUE)[,2]),asc.out$pars,d0=2),T)

asc.deg <- try(names(which(abs(asc.out$delta[,1])>=1&(asc.out$postp_d0_2[,2]>=0.99|asc.out$postp_d0_2[,3]>=0.99))),T)

}


# perform Z-test to indentify DE genes

ztest.out <- ztest(counts(cds,normalized=TRUE),conditions(cds))

ztest.deg <- rownames(ztest.out)[-log10(ztest.out$pval)>=2/abs(log2(ztest.out$FC))]


# return results

return(list(deseq.deg=deseq.deg,deseq2.deg=deseq2.deg,noiseq.deg=noiseq.deg,edger.deg=edger.deg,asc.deg=asc.deg,ztest.deg=ztest.deg))


}

