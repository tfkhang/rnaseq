##01/20/2011


############################################################
# X1 X2: two vectors of tag counts
# to estimate delta
#
#############################################################




############################################################
# Top level
#############################################################

DGE <- function(X1,X2,pars=NULL){
  a1=(20:100)^2
  a2=sapply(20:100,function(nn) sum((X1<=nn)*(X2<=nn),na.rm=T))
  xxx0=20-1+which.min(length(X1)+a1-a2)
  
  if (is.null(pars))  pars=getPars(X1,X2)
  subindex=which(X1>xxx0|X2>xxx0)

  PD1=getPostDiff(X1,X2,subindex=subindex,pars=pars)

  PDall=matrix(NA,length(X1),2)
  PDall[subindex,]=PD1$delta
  ##for small ones X1&X2<=45=xxx0

  ##kk is a interger 10^k so that X1 and X2 can be pasted without overlapping
  kk=10^(ceiling(log10(xxx0)))
  id.small=X1[-subindex]*kk+X2[-subindex]
  Xgrid=expand.grid(0:xxx0,0:xxx0)[-1,]
  id.grid=Xgrid[,1]*kk+Xgrid[,2]

  PD2 <-getPostDiff(Xgrid[,1],Xgrid[,2],pars=PD1$pars)
  PDall[-subindex,]=PD2$delta[match(id.small,id.grid)]

  rownames(PDall)=names(X1)
  list(delta=PDall,pars=PD1$pars)
}





#############################################################################################
########################            estimating prior parameters             #######################
####################################################################################

getPars <- function(x1,x2,Q1=c(0.8,0.9)){
  N1=sum(x1)
  N2=sum(x2)
  tau0=getTau(x1,x2)
  logpi=(log(pmax(1,x1)/N1)+log(pmax(1,x2)/N2))/2
  parSE=getShExp(logpi,Q1)
  pars=list("tau0"=tau0,"lambda0"=parSE[2],"rate"=parSE[1],"N1"=N1,"N2"=N2)
  pars
}

##################################################
## estimating LAMBDA AS SHIFTED EXPONENTIAL
##################################################
#1 conditional expectation of  exponential
##let x_q be the q-th quantile of exponential(a)
#E[X[X>xq]=xq+1/a
# estimating shifted exponential
getShExp <- function(x,q=c(.5,0.75)){
  xq=quantile(x,q)
  a=-diff(log(1-q))/diff(xq)#rate
  lambda0=quantile(x,q[1])+log(1-q[1])/a 
  names(a)=NULL;names(lambda0)=NULL
  c("rate"=a,"shift"=lambda0)
}

##################################################
## estimating DELTA AS NORMAL
##################################################

getTau <- function(X1,X2,q0=.995){
  N1=sum(X1);N2=sum(X2)
  lx1=log(pmax(1,X1))-log(N1)
  lx2=log(pmax(1,X2))-log(N2)
  ratio1=lx1-lx2
  xbar=(lx1+lx2)/2
  ratio1[ratio1>5]=5 #truncate it
  ratio1[-ratio1>5]=-5
  tau0=IQR(ratio1[xbar>quantile(xbar,q0)])/diff(qnorm(c(.25,.75)))
  tau0
}
#############################################################################################
########################          POSTERIOR EXPECTATION             #######################
####################################################################################
##################################################
  
getPostDiff <- function(x1,x2,subindex=NULL,pars){
  if(is.null(subindex)) subindex=1:length(x1)
  if(length(x1)!=length(x2)) stop("two libraries should have the same length")
  
  N1=pars$N1
  N2=pars$N2
  
  tau0=pars$tau0
  Ds <- seq(-15*tau0,15*tau0,0.01)
  priorD=dnorm(Ds,0,tau0)
  
  x1=x1[subindex]
  x2=x2[subindex]
  
  PostDiff <- 
    function(x1,x2,N1,N2,tau0,lambda0,alpha){
      N=N1+N2
      l1=max(log(1/N),(log(x1)+log(x2)-log(N1)-log(N2))/2)
      l2=log((x1+x2)/N)
      width1=(12*tau0-(l2-l1))*(12*tau0>(l2-l1)) # if very close, expand it
      l1=l1-width1/2;l2=l2+width1/2
      la=seq(l1,l2,len=100)+(-lambda0-l1)*(-lambda0>l1) # so now I have at least +/-6tau0 COVERAGE
#### la is the density of lambda
      
                                        #reset the range
      priorL=dexp(la+lambda0,rate=alpha)
      postD <- sapply(la,function(lamb){
        integrand=dpois(x1,N1*exp(lamb+Ds/2))*dpois(x2,N2*exp(lamb-Ds/2))*priorD
        c(ED=sum(Ds*integrand), ED2=sum(Ds^2*integrand),XgivenL=sum(integrand))})
      zzz=postD[3,]*priorL
      k=tinyIncre(zzz)
      
      
      { if(identical(k,c(1,100))){
        tmp=postD%*%priorL
        tmp1=c(tmp[1]/tmp[3],tmp[2]/tmp[3])
        tmp1[2]=sqrt(tmp1[2]-tmp1[1]^2)
        tmp1
      }
      else if(identical(k,-1)){NA}
      else{#refine it
        la=seq(la[k[1]],la[k[2]],len=100)
        priorL=dexp(la+lambda0,rate=alpha)
        k0=which(Ds>0)
        postD=sapply(la,function(lamb){
          integrand=dpois(x1,N1*exp(lamb+Ds/2))*dpois(x2,N2*exp(lamb-Ds/2))*priorD
          c(ED=sum(Ds*integrand),sum(Ds^2*integrand) ,XgivenL=sum(integrand))})
        tmp=postD%*%priorL
        tmp1=c(tmp[1]/tmp[3],tmp[2]/tmp[3])
        tmp1[2]=sqrt(tmp1[2]-tmp1[1]^2)
        tmp1
      }}
    }
  
  
  tinyIncre <- function(zzz,min.step=1e-5){# for unimodal. find the non-negligible range
    if(sum(zzz)==0){#flat
      k=-1}
    else{
      steps=diff(cumsum(zzz))/sum(zzz)
      k0=steps>min.step
      k=range(which(k0))
      k=c(max(1,k[1]),min(length(zzz),k[2]+1))
    }
    k
  }
  
##################################################  
  delta=matrix(NA,length(subindex),2)
  for(i in 1:nrow(delta)){
#cat("i=",i,"\n")
    delta[i,]=PostDiff(x1[i],x2[i],N1,N2,tau0=tau0,lambda0=-pars$lambda0,alpha=pars$rate)
    if(is.na(delta[i,1])){
      delta[i,]=c(log(x1[i]/N1)-log(pmax(1,x1[i])/N1)/2-log(pmax(1,x2[i])/N2)/2,NA)}
  }
  

  list(delta=delta,pars=pars)
}



tinyIncre <- function(zzz,min.step=1e-5){# for unimodal. find the non-negligible range
  if(sum(zzz)==0){#flat
    k=-1}
  else{
    steps=diff(cumsum(zzz))/sum(zzz)
    k0=steps>min.step
    k=range(which(k0))
    k=c(max(1,k[1]),min(length(zzz),k[2]+1))
  }
  k
}

############################################################
## call locfdr 
##  
############################################################
asc.fdr <- function(X1,X2,delta,nulltype=1){
  library(locfdr)
  N1=delta$pars$N1
  N2=delta$pars$N2
  
  pi1=X1/N1;pi2=X2/N2
  k=which((1-pi1)/(N1*pi1)<delta$pars$tau*1e-2&(1-pi2)/(N2*pi2)<delta$pars$tau*1e-2)
  D=delta$delta
  fdr=rep(NA,length(D))
  fdr[k]=locfdr(D[k],nulltype=nulltype,plot=0)$fdr
  yl=min(fdr[k][D[k]<0])
  yr=min(fdr[k][D[k]<0])
  
  fdr[-k] <- approx(D[k],fdr[k],D[-k],rule = 2)$y
  fdr
}

PostProb <- function(d,x1,x2,pars,d0=2) {
  d <- matrix(d,ncol=2)
  tmp=apply(cbind(d,x1,x2),1,function(xxx){
    callPostProb(xxx[1:2],xxx[3],xxx[4],pars,d0=d0)})
  tmp=matrix(tmp,ncol=3,byrow=T)
  colnames(tmp)=c("delta",paste("P(X1/X2>",d0,")",sep=""),
            paste("P(X2/X1>",d0,")",sep=""))
  rownames(tmp)=names(x1)
  tmp
}

callPostProb <- function(d,x1,x2,pars,d0=2){#d0 is the minimum foldchange
  if(is.na(d[2])) return(c(d[1],d[1]>0,1-d[1]>0))
  N1=pars$N1
  N2=pars$N2
  tau0=pars$tau0
  lambda0=-pars$lambda0
  alpha=pars$rate
  
  Ds <- d[1]+seq(-5*d[2],5*d[2],len=200)
  priorD=dnorm(Ds,0,tau0)
  
  N=N1+N2
  l1=max(log(1/N),(log(x1)+log(x2)-log(N1)-log(N2))/2)
  l2=log((x1+x2)/N)
  width1=(12*d[2]-(l2-l1))*(12*d[2]>(l2-l1)) # if very close, expand it
  l1=l1-width1/2;l2=l2+width1/2
  la=seq(l1,l2,len=100)+(-lambda0-l1)*(-lambda0>l1) # so now I have at least +/-6tau0 COVERAGE
#### la is the density of lambda
  
                                        #reset the range
  priorL=dexp(la+lambda0,rate=alpha)
  k1=which(Ds>log(d0))
  k2=which(-Ds>log(d0))
  postD=sapply(la,function(lamb){
    integrand=dpois(x1,N1*exp(lamb+Ds/2))*dpois(x2,N2*exp(lamb-Ds/2))*priorD
    c(ED=sum(Ds*integrand),p1=sum(integrand[k1]),p2=sum(integrand[k2]),XgivenL=sum(integrand))})
  zzz=postD[4,]*priorL
  k=tinyIncre(zzz)

  { if(identical(k,c(1,100))){
    tmp=postD%*%priorL
    tmp[1:3]/tmp[4]

  }
  else if(identical(k,-1)){NA}
  else{#refine it
    la=seq(la[k[1]],la[k[2]],len=100)
    priorL=dexp(la+lambda0,rate=alpha)
    postD=sapply(la,function(lamb){
      integrand=dpois(x1,N1*exp(lamb+Ds/2))*dpois(x2,N2*exp(lamb-Ds/2))*priorD
      c(ED=sum(Ds*integrand),p1=sum(integrand[k1]),p2=sum(integrand[k2]),XgivenL=sum(integrand))})
    
    tmp=postD%*%priorL
    tmp[1:3]/tmp[4]
    
  }}
}
