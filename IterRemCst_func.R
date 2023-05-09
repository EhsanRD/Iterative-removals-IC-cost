#setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Iterative-Removals-SW\\IterativeInfContent")
source("CRTVarAdj_func.R", local=TRUE)
library("shiny")
library("ggplot2")
library("reshape2")
library("plyr")
library("swCRTdesign")
library("matrixcalc")
library("scales")
library("tidyverse")
library("shinythemes")
library("Matrix")
library("plotly")
library("RColorBrewer")
library("tidyr")
#generating design matrices.
SWdesmat <- function(Tp) {
  Xsw <- matrix(data=0, ncol = Tp, nrow = (Tp-1))
  for(i in 1:(Tp-1)) {
    Xsw[i,(i+1):Tp] <- 1
  }
  return(Xsw)
}

DsgnCst <- function(Xdes,Tp, N,m, c, s,sprim, R, Rprim){
  # Returns total trial cost, given:
  #   - number of periods, Tp
  #   - number of clusters per each sequence, N
  #   - number of subjects per cluster-period, m
  #   - cost per cluster, c
  #   - cost per subject under intervention condition, s
  #   - cost per subject under intervention condition, sprim
  #   - restart cost under intervention condition, R
  #   - restart cost under intervention condition, Rprim
  #   - number of clusters with at least one non-missing cell, K
  #   - total non-missing cells, o
  #   - length of gaps, l
  #   - number of gaps, n
 
  #assume one cluster per each sequence
  n_1=numeric(Tp-1)
  len_gaps_1=numeric(Tp-1)
  n_2=numeric(Tp-1)
  len_gaps_2=numeric(Tp-1)
  k=numeric(Tp-1)
  #define custom function to calculate min
  custom_min <- function(x) {if (length(x)>0) min(x) else Inf}
  custom_max <- function(x) {if (length(x)>0) max(x) else Inf}
  
  
  for (i in 1:(Tp-1)){
    #length of gaps starting and ending with intervention condition
    l_1 <- custom_min(which((Xdes[i,]==1))) < which(is.na(Xdes[i,])) &
      which(is.na(Xdes[i,]))  < custom_max(which((Xdes[i,]==1)))
    #length of gaps starting and ending with control condition
    l_2 <- custom_min(which((Xdes[i,]==0))) < which(is.na(Xdes[i,])) &
      which(is.na(Xdes[i,]))  < custom_max(which((Xdes[i,]==0)))
    
    len_gaps_1[i]=sum(l_1)
    len_gaps_2[i]=sum(l_2)
    #identify the number of gaps
    n_1[i]=ifelse(sum(l_1)>0,1,0)
    n_1[i]=ifelse(sum(diff(which(is.na(Xdes[i,]))))!=sum(l_1)-1 & sum(l_1)>1,sum(diff(which(is.na(Xdes[i,]))[l_1==TRUE])>1)+1,n_1[i])
    
    n_2[i]=ifelse(sum(l_2)>0,1,0)
    n_2[i]=ifelse(sum(diff(which(is.na(Xdes[i,]))))!=sum(l_2)-1 & sum(l_2)>1,sum(diff(which(is.na(Xdes[i,]))[l_2==TRUE])>1)+1,n_2[i])
    #count clusters that have at least one non-missing value.
    k[i]=ifelse(sum(!is.na(Xdes[i,]))>0,1,0)
  }
  #consider equal number of clusters per each sequence
  o=sum(!is.na(Xdes[,]))*N

  R_1=sum(n_1)*R*N

  R_2=sum(n_2)*Rprim*N

  k=sum(k)

  cvec <- k*N*c + 0.5*o*m*(s+sprim) + R_1 + R_2
  
  return(cvec)
}

#assume one cluster is randomised to each sequence. 
#function to calculate the information content considering centrosymmetric property. 
ICPair = function(Xdes,m,rho0,r,type) {
  
  Tp <- ncol(Xdes)
  K  <- nrow(Xdes)
  varD=CRTVarGeneralAdj(Xdes,m,rho0,r,type)
  ICmat<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
  
  for(i in 1:nrow(Xdes)){
    for (j in 1:ncol(Xdes)){
      if(is.na(Xdes[i,j])==TRUE | is.na(Xdes[K-i+1,Tp-j+1])==TRUE){
        ICmat[i,j] <- NA
        ICmat[K-i+1,Tp-j+1] <- NA
      }
      else if(is.na(Xdes[i,j])==FALSE & is.na(Xdes[K-i+1,Tp-j+1])==FALSE) {
        Xdesij <- Xdes
        Xdesij[i,j] <- NA
        Xdesij[K-i+1,Tp-j+1] <- NA
        ICmat[i,j] <- round(CRTVarGeneralAdj(Xdesij,m,rho0,r,type)/varD,10)
        ICmat[K-i+1,Tp-j+1] <- ICmat[i,j]
        
        if (is.na(ICmat[i,j])==TRUE & is.na(ICmat[K-i+1,Tp-j+1]==TRUE)) {
          ICmat[i,j] <-101.101
          ICmat[K-i+1,Tp-j+1] <- ICmat[i,j]
        }
      }
      #avoid replicating loop
      if (i*j==(nrow(Xdes)*ncol(Xdes)/2)){
        break
      }
    }
  }
  return(ICmat)
}

IterRemove = function(Tp,N,m,rho0,r,type,c,s,sprim,R,Rprim){
  
  K=Tp-1
  mval <- list()    #minimum lowest information content values
  Xdlist <- list()  #design matrix
  dlist <- list()   #information content matrix
  Xdlist[[1]] <- SWdesmat(Tp)
  dlist[[1]] <- ICPair(Xdlist[[1]],m,rho0,r,type)
  cvec<- c()
  cvec[1]<-DsgnCst(Xdlist[[1]],Tp,N, m, c, s,sprim, R,Rprim)[[1]]
  varvec<- c()
  varvec[1]<-CRTVarGeneralAdj(Xdlist[[1]],m,rho0,r,type)
  #removal of pairs of cells
  for (i in 2:(Tp*K/2-1)){ #most minimal design 
    mval <- which(dlist[[i-1]]==min(dlist[[i-1]],na.rm = TRUE), arr.ind = TRUE)
    Xdlist[[i]]=Xdlist[[i-1]]
    #remove a pair of centrosymmetric cells
    #re-order indices by cluster and period
    mval <- mval[order(mval[,1],mval[,2]),]
    cellid <- mval[1,]
    clustid<-cellid[1]
    perid<-cellid[2]
    Xdlist[[i]][clustid,perid]<- NA
    Xdlist[[i]][K+1-clustid,Tp+1-perid]<- NA
    
    cvec[i]=DsgnCst(Xdlist[[i]],Tp,N, m, c, s,sprim, R,Rprim)
    varvec[i]<-CRTVarGeneralAdj(Xdlist[[i]],m,rho0,r,type)
    
    dlist[[i]] = ICPair(Xdlist[[i]],m,rho0,r,type)
    
  }
  return(list(dlist, Xdlist, varvec, cvec))
}



