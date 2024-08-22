pow <- function(vars, effsiz, siglevel=0.05){
  z <- qnorm(siglevel/2)
  pow <- pnorm(z + sqrt(1/vars)*effsiz)
  return(pow)
}

DsgnCst <- function(Xdes,Tp, N,m, c, p,pprim, g, gprim){
  # Returns total trial cost, given:
  #   - number of periods, Tp
  #   - number of clusters per each sequence, N
  #   - number of subjects per cluster-period, m
  #   - cost per cluster, c
  #   - cost per subject under intervention condition, p
  #   - cost per subject under control condition, pprim
  #   - restart cost under intervention condition, g
  #   - restart cost under control condition, gprim
  #   - number of sequences with at least one non-missing cell, S
  #   - total non-missing cells, o
  #   - length of gaps, l_1, l_2
  #   - number of gaps, n_1, n_2
  
  #assume one cluster per each sequence
  n_1=numeric(Tp-1)
  len_gaps_1=numeric(Tp-1)
  n_2=numeric(Tp-1)
  len_gaps_2=numeric(Tp-1)
  Is=numeric(Tp-1)
  Ts_1=numeric(Tp-1)
  Ts_2=numeric(Tp-1)
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
    #count sequences that have at least one non-missing value.
    Is[i]=ifelse(sum(!is.na(Xdes[i,]))>0,1,0)
    Ts_1[i]=sum(Xdes[i,] == 1, na.rm = TRUE)
    Ts_2[i]=sum(Xdes[i,] == 0, na.rm = TRUE)
  }
  #consider equal number of clusters per each sequence
  Sum_Is=sum(Is)
  Sum_Ts_1=sum(Ts_1)
  Sum_Ts_2=sum(Ts_2)
  Sum_n_1=sum(n_1)
  Sum_n_2=sum(n_2)
  
  cvec <- c*Sum_Is + m*p*Sum_Ts_1 + m*pprim*Sum_Ts_2 + g*Sum_n_1 + gprim*Sum_n_2
  cvec <- N*cvec
  
  return(cvec)
}


#Calculate the cost efficiency metric
CEcal = function(Xdes,N,m,rho0,r,type,c,p,pprim,g,gprim)  {
  
  Tp <- ncol(Xdes)
  S  <- nrow(Xdes)
  varD = CRTVarGeneralAdj(Xdes,m,rho0,r,type)/N
  #cvecD = DsgnCst(Xdes,Tp,N, m, c, p,pprim, g,gprim)
 
  CEmat<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
  #REmat<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
  
  for(i in 1:nrow(Xdes)){
    for (j in 1:ncol(Xdes)){
      if(is.na(Xdes[i,j])==TRUE){
        CEmat[i,j] <- NA
        #REmat[i,j] <- NA
      }
      else if(is.na(Xdes[i,j])==FALSE) {
        Xdesij <- Xdes
        Xdesij[i,j] <- NA
        #N will be cancel out for equal number of clusters per sequences
        CEmat[i,j] <- round((1/CRTVarGeneralAdj(Xdesij,m,rho0,r,type)/N)/
                              DsgnCst(Xdesij,Tp, N,m, c, p,pprim, g, gprim),10)
        #REmat[i,j] <-  CEmat[i,j] / cvecD
        if (is.na(CEmat[i,j])==TRUE) {
          CEmat[i,j] <-101.101
        }
      }
      #avoid replicating loop
      # if (i*j==(nrow(Xdes)*ncol(Xdes)/2)){
      #   break
      # }
    }
  }
  return(CEmat)
  #return(list(CEmat,REmat))
}

#Remove the design that has the highest CE with removing the cell

IterRemove_CE = function(Tp,N,m,rho0,r,type,c,p,pprim,g,gprim,effsiz,accept_pwr){
  
  S=Tp-1
  mval <- list()    #minimum lowest CE
  Xdlist <- list()  #design matrix
  dlist <- list()   #cost efficiency matrix
  Xdlist[[1]] <- SWdesmat(Tp)
  dlist[[1]] <- CEcal(Xdlist[[1]],N,m,rho0,r,type,c,p,pprim,g,gprim)
 
  cvec<- numeric()
  cvec[1]<-DsgnCst(Xdlist[[1]],Tp, N,m, c, p,pprim, g, gprim)[[1]]
  varvec<- numeric()
  varvec[1]<-CRTVarGeneralAdj(Xdlist[[1]],m,rho0,r,type)/N
  pwvec<- numeric()
  pwvec[1]<- pow(varvec[1],effsiz,siglevel=0.05)*100
  if(accept_pwr < pwvec[1]) {
    #removal of single cells
    for (i in 2:(Tp*S-2)){ #most minimal design
      mval <- which(dlist[[i-1]]==max(dlist[[i-1]],na.rm = TRUE), arr.ind = TRUE)
      Xdlist[[i]]=Xdlist[[i-1]]
      #remove a single cell
      #re-order indices by cluster and period
      #mval <- mval[order(mval[,1],mval[,2]),]
      #force R to return a matrix even when it has only one row
      mval <- mval[order(mval[,1], mval[,2]), , drop = FALSE]
      cellid <- mval[1,]
      clustid<-cellid[1]
      perid<-cellid[2]
      Xdlist[[i]][clustid,perid]<- NA
      #Xdlist[[i]][S+1-clustid,Tp+1-perid]<- NA
      cvec[i]=DsgnCst(Xdlist[[i]],Tp,N, m, c, p,pprim, g,gprim)
      varvec[i]<-CRTVarGeneralAdj(Xdlist[[i]],m,rho0,r,type)/N
      pwvec[i]<- pow(varvec[i],effsiz,siglevel=0.05)*100
      
      if (accept_pwr < pwvec[i]) {
        dlist[[i]] = CEcal(Xdlist[[i]],N,m,rho0,r,type,c,p,pprim,g,gprim)
      }
      else {
        Xdlist <- Xdlist[-i]
        cvec <- cvec[-i]
        varvec <- varvec[-i]
        pwvec <- pwvec[-i]
        break
      }
      
    }
    return(list(dlist, Xdlist, varvec, cvec,pwvec))
  }
  else {
    return(NULL)
  }
}

