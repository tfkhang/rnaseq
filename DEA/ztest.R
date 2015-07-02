
ztest <- function(data.table,cond){

whichx <- names(cond)[cond==levels(cond)[1]]
whichy <- names(cond)[cond==levels(cond)[2]]
nx <- sum(data.table[,whichx])
ny <- sum(data.table[,whichy])

as.data.frame(t(apply(data.table,1,function(f){
x <- sum(f[whichx])
y <- sum(f[whichy])
FC <- y/x
test.stat <- sum(x,y)/sum(nx,ny)
test.stat <- (x/nx-y/ny)/sqrt(test.stat*(1-test.stat)/sum(nx,ny))
pval <- pnorm(-abs(test.stat))*2
return(c(FC=FC,pval=pval))})))

}
