# Functions for generating design matrices
SWdesmat <- function(Tp) {
  Xsw <- matrix(data=0, ncol = Tp, nrow = (Tp-1))
  for(i in 1:(Tp-1)) {
    Xsw[i,(i+1):Tp] <- 1
  }
  return(Xsw)
}
##############
CRTVarGeneralAdj <- function(Xmat, m, rho0, r, type) {
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  
    Tp <- ncol(Xmat)
    K  <- nrow(Xmat)
    Xvec <- as.vector(t(Xmat))
    
    stackI <- matrix(rep(diag(1,Tp)), nrow=K*Tp, ncol=Tp, byrow=TRUE)
    Zmat <- cbind(stackI, Xvec)
    Zmat <- Zmat[!is.na(Xvec),]
    
    dropind=which(colSums(!is.na(Xmat)) == 0)
    
    if(length(dropind>0)){
      Zmat <- Zmat[,-dropind]
    }
    #Constant decay
    if(type==0) { 
      Vi <-diag(sig2 +(1-r)*sig2CP, Tp) + matrix(data=sig2CP*r, nrow=Tp, ncol=Tp)
    }
    #exponential decay structure
    if(type==1) {
      Vi <- diag(sig2,Tp) + sig2CP*(r^abs(matrix(1:Tp,nrow=Tp, ncol=Tp, byrow=FALSE) -
                          matrix(1:Tp,nrow=Tp, ncol=Tp, byrow=TRUE)))
    }
  #Variance matrix for all clusters
  Vall <- kronecker(diag(1,K), Vi)
  Vall <- Vall[!is.na(Xvec),!is.na(Xvec)]
  
  #there will be problems if Zmat is not of full column rank
  if(rankMatrix(Zmat)[1] < ncol(Zmat)) return(NA)
  else return(solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)])
}

