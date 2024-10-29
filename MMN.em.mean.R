library(mvtnorm)
#library(clusterGeneration)

mmn.em_mean<-function(y,k,memb=NULL,q=4,ort=TRUE,test=NULL,it=1000,eps=0.0000001,seme=7,val.init=NULL,init="momenti",model.name=c("VVV","VVV","U"),
                      m=1)
{
  
  ## model.name can be:
  ## [1] PHI: - temporal: GAR, GARI (GAR+isotropic for D), EGAR, EGARI, with m=1,...,T-1
  ##          - or not: VVV,EEV, EEE, III, VVI, EII, EEI
  ## [2] OMEGA: VVV,EEV, EEE, III, VVI, EII, EEI  (Vedere paper su STAT&Comp)
  ## y ? di dimensione r x p x numobs per ottenere (r sarebbe T)
  
  
  ptm = proc.time()
  set.seed(seme)
  if(!is.null(test)){
    y.test <- y[,,test]
    y <- y[,,-test]
  } 
  r=dim(y)[1]
  p=dim(y)[2]
  n=dim(y)[3]
  q <- min(q, r-1)
  
  yy=t(apply(y,3,c))
  
  if (!is.null(val.init)) {M=val.init$M
  U=val.init$U
  V=val.init$V
  CC=V
  A=U
  w=val.init$w
  lambda=w
  csi=w
  }
  
  
  ######################################
  ###         INIZIALIZZAZIONE       ###
  ######################################
  
  
  ## inizializzazione con i momenti di popolazione
  if (is.null(val.init)) {
    
    M=array(0,c(r,p,k))
    mu=matrix(0,r,k)
    nu=matrix(0,p,k)
    alfa=matrix(0,r,k)
    beta=matrix(0,p,k)
    delta=rep(0,k)
    if(model.name[3] == "PI" | model.name[3] == "P" | model.name[3] == "PIO"){
      Q <- cbind(1/sqrt(r),poly(1:r,q-1,raw = !ort,simple=TRUE))
    }
    a=matrix(0,q,k)
    b=matrix(0,q,k)
    P=array(0,dim =c(q,p,k))
    
    if (k>1) {
      if(is.null(memb)) memb<-kmeans(yy,k)$cl
      #if(n>r*p) {memb <- summary(Mclust(yy, modelName = "VVV", G = k))$classification}
      #else {
      #memb <- summary(Mclust(yy, modelName = "VVI", G = k))$classification#}
    }
    else {memb=rep(1,n)}
    for  (i in 1:k) {
      if ((table(memb)[i])<2) {
        memb[sample(1:n,2,replace=FALSE)]=i;  cat("table(memb)[i]<2 \n")
      }###check per componenti con 1 solo elemento
    }
    
    for (i in 1:k){
      M[,,i]=apply(y[,,memb==i],c(1,2),mean)
      if(model.name[3] == "A"){
        nu[,i]  <- apply(y[,,memb==i], 2, mean)
        mu[,i] <- apply(y[,,memb==i], 1, mean) + nu[1,i]
        nu[,i] <- nu[,i] - nu[1,i]
        M[,,i] <- matrix(mu[,i], r, p) + matrix(nu[,i], r, p, byrow = TRUE) 
      }
      
      if(model.name[3] == "I"){
        nu[,i]  <- apply(y[,,memb==i], 2, mean)
        mu[,i] <- apply(y[,,memb==i], 1, mean)
        delta[i] <- nu[1,i] + mu[1,i]
        nu[,i] <- nu[,i] - nu[1,i]
        mu[,i] <- mu[,i] - mu[1,i]
        #M[,,i]  <- matrix(mu[,i], r, p) + matrix(nu[,i], r, p, byrow = TRUE) + delta
        M[,,i] <- matrix(mu[,i], r, p) + matrix(nu[,i], r, p, byrow = TRUE) + mu[,i]%*%t(nu[,i]) + delta[i]
      }
      if(model.name[3] == "GI"){
        delta[i] <- M[1,1,i]
        nu[,i] <- M[1,,i] - delta[i]
        mu[,i] <- M[,1,i] - delta[i]
        SVD <- svd(M[,,i] - ( matrix(mu[,i], r, p) + matrix(nu[,i], r, p, byrow = TRUE) + delta[i]))  
        alfa[,i] <- SVD$u[,1]*sqrt(SVD$d[1])
        alfa[,i] <- alfa[,i] - alfa[1,i]
        alfa[,i] <- alfa[,i]/alfa[2,i]
        beta[,i] <- SVD$v[,1]*sqrt(SVD$d[1])
        beta[,i] <- beta[,i] - beta[1,i]
        beta[,i] <- beta[,i]/beta[2,i]
        #M[,,i] <- M[,,i] + alfa[,i] %*% t(beta[,i])
        M[,,i] <- matrix(mu[,i], r, p) + matrix(nu[,i], r, p, byrow = TRUE) + alfa[,i]%*%t(beta[,i]) + delta[i]
      }
      if(model.name[3] == "PI"){
        delta[i] <- M[1,1,i]
        nu[,i] <- M[1,,i] - delta[i]
        mu[,i] <- M[,1,i] - delta[i]
        SVD <- svd(M[,,i] - ( matrix(mu[,i], r, p) + matrix(nu[,i], r, p, byrow = TRUE) + delta[i])) 
        alfa[,i] <- SVD$u[,1]*sqrt(SVD$d[1])
        alfa[,i] <- alfa[,i] - alfa[1,i]
        alfa[,i] <- alfa[,i]/alfa[2,i]
        beta[,i] <- SVD$v[,1]*sqrt(SVD$d[1])
        beta[,i] <- beta[,i] - beta[1,i]
        beta[,i] <- beta[,i]/beta[2,i]
        #M[,,i] <- M[,,i] + alfa[,i] %*% t(beta[,i])
        b[,i] <- lm(alfa[,i] ~ Q - 1)$coefficients
        a[,i] <- lm(mu[,i] ~ Q - 1)$coefficients
        alfa[-(1:2),i] <- Q[-(1:2),]%*%b[,i]
        M[,,i] <- matrix(c(0,Q[-1,]%*%a[,i]), r, p) + matrix(nu[,i], r, p, byrow = TRUE) + c(0,1,Q[-(1:2),]%*%b[,i])%*%t(beta[,i]) + delta[i]
      }
      if(model.name[3] == "PIO"){
        nu[,i]  <- apply(y[,,memb==i], 2, mean)
        nu[,i] <- nu[,i] - nu[1,i]
        SVD <- svd(M[,,i] - matrix(nu[,i], r, p, byrow = TRUE)) 
        beta[,i] <- SVD$v[,1]*sqrt(SVD$d[1])
        alfa[,i] <- SVD$u[,1]*sqrt(SVD$d[1])*beta[1,i]
        beta[,i] <- beta[,i]/beta[1,i]
        #M[,,i] <- M[,,i] + alfa[,i] %*% t(beta[,i])
        b[,i] <- lm(alfa[,i] ~ Q - 1)$coefficients
      }
    }
    if (init=="rand") for (i in 1:k) {M[,,i]=M[,,i]+matrix(runif(r*p,-10,10),r,p)}
    
    U=array(0,c(r,r,k))
    V=array(0,c(p,p,k))
    
    if (model.name[1]=="VVV" | model.name[1]=="EEV") for (i in 1:k) U[,,i]=var(scale(apply(y[,,memb==i],c(1),c)))+genPositiveDefMat(r,covMethod="unifcorrmat",rangeVar=c(0,2))$Sigma
    if (model.name[2]=="VVV" | model.name[2]=="EEV") for (i in 1:k) V[,,i]=var(apply(y[,,memb==i],c(2),c))+genPositiveDefMat(p,covMethod="unifcorrmat",rangeVar=c(0,2))$Sigma
    
    Var1=var(scale(apply(y,c(1),c)))+genPositiveDefMat(r,covMethod="unifcorrmat",rangeVar=c(0,1))$Sigma
    Var2=var(apply(y,c(2),c))+genPositiveDefMat(p,covMethod="unifcorrmat",rangeVar=c(0,1))$Sigma
    
    if (model.name[1]=="EEE") for (i in 1:k) U[,,i]=Var1
    if (model.name[2]=="EEE") for (i in 1:k) V[,,i]=Var2
    
    if (model.name[1]=="III") for (i in 1:k) U[,,i]=diag(r)
    if (model.name[2]=="III") for (i in 1:k) V[,,i]=diag(p)
    
    lambda=matrix(0,r)
    csi=matrix(0,p)
    rho=matrix(0,r)
    A=array(0,c(r,r,k))
    D=array(0,c(r,r,k))
    T=array(0,c(r,r,k))
    CC=array(0,c(p,p,k))
    
    if (model.name[1]=="VVI") for (i in 1:k) {U[,,i]=var(scale(apply(y[,,memb==i],c(1),c)))+genPositiveDefMat(r,covMethod="unifcorrmat",rangeVar=c(0,1))$Sigma
    lambda[i]=det(U[,,i])^(1/r)
    A[,,i]=U[,,i]/lambda[i]
    }
    
    if (model.name[2]=="VVI") for (i in 1:k) {V[,,i]=var(apply(y[,,memb==i],c(2),c))+genPositiveDefMat(p,covMethod="unifcorrmat",rangeVar=c(0,1))$Sigma
    csi[i]=det(V[,,i])^(1/p)
    CC[,,i]=V[,,i]/csi[i]}
    
    if (model.name[1]=="EEI") for (i in 1:k) {U[,,i]=Var1
    lambda[i]=det(U[,,i])^(1/r)
    A[,,i]=U[,,i]/lambda[i]
    }
    
    if (model.name[2]=="EEI") for (i in 1:k) {V[,,i]=Var2
    csi[i]=det(V[,,i])^(1/p)
    CC[,,i]=V[,,i]/csi[i]
    }
    
    
    if (model.name[1]=="VII") for (i in 1:k) {U[,,i]=var(scale(apply(y[,,memb==i],c(1),c)))+genPositiveDefMat(r,covMethod="unifcorrmat",rangeVar=c(0,1))$Sigma
    lambda[i]=det(U[,,i])^(1/r)
    A[,,i]=diag(r)
    U[,,i]=lambda[i]*A[,,i]
    }
    
    if (model.name[2]=="VII") for (i in 1:k) {V[,,i]=var(apply(y[,,memb==i],c(2),c))+genPositiveDefMat(p,covMethod="unifcorrmat",rangeVar=c(0,1))$Sigma
    csi[i]=det(V[,,i])^(1/p)
    CC[,,i]=diag(p)
    V[,,i]=csi[i]*CC[,,i]
    }
    
    
    
    if (model.name[1]=="EII") for (i in 1:k) {U[,,i]=Var1
    lambda[i]=det(U[,,i])^(1/r)
    A[,,i]=diag(r)
    U[,,i]=lambda[i]*A[,,i]
    }
    
    if (model.name[2]=="EII") for (i in 1:k) {V[,,i]=Var2
    csi[i]=det(V[,,i])^(1/p)
    CC[,,i]=diag(p)
    V[,,i]=csi[i]*CC[,,i]
    }
    
    if (model.name[1]=="GAR") for (i in 1:k) {  T[,,i]=matrix(runif(r*r,-1,0),r,r)
    T[,,i][upper.tri(T[,,i],diag=TRUE)]<-0
    for (h in 1:r) if ((h+m+1)<=r) T[(h+m+1):r,h,i]<-0
    diag(T[,,i])=1
    D[,,i]<-diag(runif(r))
    U[,,i]=t(T[,,i])%*%ginv(D[,,i])%*%T[,,i]
    U[,,i]=ginv(U[,,i])
    }
    
    if (model.name[1]=="GARI") for (i in 1:k) { T[,,i]=matrix(runif(r*r,-1,0),r,r)
    T[,,i][upper.tri(T[,,i],diag=TRUE)]<-0
    for (h in 1:r) if ((h+m+1)<=r) T[(h+m+1):r,h,i]<-0
    diag(T[,,i])=1
    D[,,i]<-diag(rep(runif(1),r))
    U[,,i]=t(T[,,i])%*%ginv(D[,,i])%*%T[,,i]
    U[,,i]=ginv(U[,,i])
    }
    
    
    if (model.name[1]=="EGAR")  {
      TT<-matrix(runif(r*r,-1,0),r,r)
      DD<-diag(runif(r))
      for (i in 1:k) {
        T[,,i]=TT
        T[,,i][upper.tri(T[,,i],diag=TRUE)]<-0
        for (h in 1:r) if ((h+m+1)<=r) T[(h+m+1):r,h,i]<-0
        diag(T[,,i])=1
        D[,,i]<-DD
        U[,,i]=t(T[,,i])%*%ginv(D[,,i])%*%T[,,i]
        U[,,i]=ginv(U[,,i])
      }
    }
    
    
    
    if (model.name[1]=="EGARI")  {
      TT<-matrix(runif(r*r,-1,0),r,r)
      DD<-diag(rep(runif(1),r))
      for (i in 1:k) {
        T[,,i]=TT
        T[,,i][upper.tri(T[,,i],diag=TRUE)]<-0
        for (h in 1:r) if ((h+m+1)<=r) T[(h+m+1):r,h,i]<-0
        diag(T[,,i])=1
        D[,,i]<-DD
        U[,,i]=t(T[,,i])%*%ginv(D[,,i])%*%T[,,i]
        U[,,i]=ginv(U[,,i])
      }
    }
    
    
    
    
    ## inizializzazione dei pesi casuale
    w=table(memb)/n
    
    if (init=="rand") {w=runif(k)
    w=w/sum(w)
    w=t(w)}
    
  }
  
  
  ######################################
  ###             INIZIO EM          ###
  ######################################
  
  
  f.y.z=matrix(0,n,k)
  for (i in 1:k) f.y.z[,i]=dmm(y,M[,,i],U[,,i],V[,,i], log = TRUE)
  f.y.z=ifelse(is.na(f.y.z),mean(f.y.z, na.rm=T),f.y.z)
  f.z.y=matrix(1,n,k)
  
  likelihood=NULL
  hh=0
  ratio=1000
  lik=-.Machine$double.xmax
  
  B=matrix(0,r,r)
  DD=array(0,c(r,r,k))
  W=matrix(0,p,p)
  L=array(0,c(p,p,k))
  omega1=matrix(0,r,k)
  omega2=matrix(0,p,k)
  
  while (hh < 10 |((hh < it) & (ratio > eps ))) {
    hh<-hh+1
    
    w.old=w
    V.old=V
    U.old=U
    M.old=M
    f.z.y.old=f.z.y
    
    ######## E-STEP
    
    medie<-colMeans(f.y.z)
    f.y.z.star = f.y.z
    fw <- matrix(w,n,k,byrow=TRUE)*exp(f.y.z)
    
    if(k > 1){
      
      if(sum(rowSums(fw)==0)>0){
        #cat(paste0("iter ", hh, ": sum(rowSums(fw)==0)>0 \n"))
        f.y.z.star[which(rowSums(fw)==0),] = 
          f.y.z.star[which(rowSums(fw)==0),] -
          apply(matrix(log(w),sum(rowSums(fw)==0),k,byrow=TRUE) + 
                  (f.y.z.star[which(rowSums(fw)==0),]),
                1,min) +
          log(.Machine$double.xmin/
                as.numeric(w[apply(
                  matrix(log(w),sum(rowSums(fw)==0),k,byrow=TRUE)
                  +(f.y.z.star[which(rowSums(fw)==0),]),
                  1,which.min)]))
      }
      
      if(sum(f.y.z.star > log(.Machine$double.xmax))>0){
        #cat(paste0("iter ", hh, ": sum(f.y.z.star > log(.Machine$double.xmax))>0 \n"))
        f.y.z.star[which(rowSums(f.y.z.star > log(.Machine$double.xmax))>0),] = 
          f.y.z.star[which(rowSums(f.y.z.star > log(.Machine$double.xmax))>0),] -
          f.y.z.star[which(rowSums(f.y.z.star > log(.Machine$double.xmax))>0),] + 
          log(.Machine$double.xmax)
      }
      
      f.y.z.star = exp(f.y.z.star)
      
      for (i in 1:k) f.z.y[,i]=f.y.z.star[,i]*w[i]/(w%*%t(f.y.z.star))
      
      f.z.y=ifelse(is.na(f.z.y),mean(f.z.y, na.rm=T),f.z.y)
    }
    
    ######## M-STEP
    
    for (i in 1:k) {
      if(model.name[3]=="U"){
        num0=0
        for(j in 1:n) num0=num0+f.z.y[j,i]*y[,,j]
        M[,,i]=num0/sum(f.z.y[,i])
      }
      
      if(model.name[3]=="A"){
        nu[,i] <- colSums(matrix(f.z.y[,i],n,p)*apply(y - array(mu[,i], dim = c(r, p, n)), c(3,2), mean))/sum(f.z.y[,i])
        nu[,i] <- nu[,i] - nu[1,i]
        mu[,i] <- colSums(matrix(f.z.y[,i],n,r)*apply(y - array(matrix(nu[,i], r, p, byrow = TRUE), dim = c(r,p,n)), c(3,1), mean))/sum(f.z.y[,i])
        M[,,i] <- matrix(mu[,i], r, p) + matrix(nu[,i], r, p, byrow = TRUE)
      }
      
      if(model.name[3]=="I"){
        nu[,i] <- 0
        for(j in 1:n){
          nu[-1,i] <- nu[-1,i]+
            f.z.y[j,i]* 
            ginv(ginv(V[,,i])[-1,-1])%*%
            ginv(V[,,i])[-1,]%*%
            (t(y[,,j]) - matrix(mu[,i], p, r, byrow = TRUE) - delta[i])%*%
            ginv(U[,,i])%*%
            (1+mu[,i])
        }
        
        nu[-1,i] <- nu[-1,i]/sum(f.z.y[,i])/((t(1+mu[,i])%*%ginv(U[,,i])%*%(1+mu[,i]))[1,1])
        
        mu[,i] <- 0
        for(j in 1:n){
          mu[-1,i] <- mu[-1,i]+
            f.z.y[j,i]* 
            ginv(ginv(U[,,i])[-1,-1])%*%
            ginv(U[,,i])[-1,]%*%
            (y[,,j] - matrix(nu[,i], r, p, byrow = TRUE) - delta[i])%*%
            ginv(V[,,i])%*%
            (1+nu[,i])
        }
        mu[-1,i] <- mu[-1,i]/sum(f.z.y[,i])/((t(1+nu[,i])%*%ginv(V[,,i])%*%(1+nu[,i]))[1,1])
        
        delta[i] <- 0
        for(j in 1:n){
          delta[i] <- delta[i] + 
            f.z.y[j,i]*
            sum(
              ginv(V[,,i])%*%
                (t(y[,,j]) - matrix(nu[,i], p, r) - matrix(mu[,i], p, r, byrow = TRUE) - nu[,i]%*%t(mu[,i]))%*%
                ginv(U[,,i])
            )
        }
        delta[i] <- delta[i]/sum(f.z.y[,i])/sum(ginv(U[,,i]))/sum(ginv(V[,,i]))
        
        M[,,i] <- matrix(mu[,i], r, p) + matrix(nu[,i], r, p, byrow = TRUE) + mu[,i]%*%t(nu[,i]) + delta[i]
      }
      
      if(model.name[3]=="GI"){
        nu[,i] <- 0
        
        for(j in 1:n){
          nu[-1,i] <- nu[-1,i]+
            f.z.y[j,i]* 
            ginv(V[,,i])[-1,]%*%
            (t(y[,,j]) - matrix(mu[,i], p, r, byrow = TRUE) - beta[,i]%*%t(alfa[,i]) - delta[i]
            )%*%
            ginv(U[,,i])%*%
            (rep(1,r))
        }
        
        nu[-1,i] <- ginv(ginv(V[,,i])[-1,-1])%*%nu[-1,i]/sum(f.z.y[,i])/sum(ginv(U[,,i]))
        
        mu[,i] <- 0
        
        for(j in 1:n){
          mu[-1,i] <- mu[-1,i]+
            f.z.y[j,i]* 
            ginv(U[,,i])[-1,]%*%
            (y[,,j] - matrix(nu[,i], r, p, byrow = TRUE) - alfa[,i]%*%t(beta[,i]) - delta[i])%*%
            ginv(V[,,i])%*%
            (rep(1,p))
        }
        mu[-1,i] <- ginv(ginv(U[,,i])[-1,-1])%*%mu[-1,i]/sum(f.z.y[,i])/sum(ginv(V[,,i]))
        
        beta[,i] <- c(0,1,rep(0,p-2))
        for(j in 1:n){
          beta[-(1:2),i] <- beta[-(1:2),i] + 
            f.z.y[j,i]*ginv(V[,,i])[-(1:2),]%*%
            (t(y[,,j]) - matrix(mu[,i], p, r, byrow = TRUE) - matrix(nu[,i], p, r) - delta[i]
            )%*%
            ginv(U[,,i])%*%
            alfa[,i] 
        }
        beta[-(1:2),i] <-  ginv(ginv(V[,,i])[-(1:2),-(1:2)])%*%
          (beta[-(1:2),i]/sum(f.z.y[,i])/((t(alfa[,i])%*%ginv(U[,,i])%*%alfa[,i])[1,1])-ginv(V[,,i])[-(1:2),2])
        
        alfa[,i] <- c(0,1,rep(0,r-2))
        for(j in 1:n){
          alfa[-(1:2),i] <- alfa[-(1:2),i] + 
            f.z.y[j,i]* ginv(U[,,i])[-(1:2),]%*%
            (y[,,j] - matrix(mu[,i], r, p) - matrix(nu[,i], r, p, byrow = TRUE) - delta[i]
            )%*%
            ginv(V[,,i])%*%
            beta[,i]
        }
        alfa[-(1:2),i] <-  ginv(ginv(U[,,i])[-(1:2),-(1:2)])%*%
          (alfa[-(1:2),i]/sum(f.z.y[,i])/((t(beta[,i])%*%ginv(V[,,i])%*%beta[,i])[1,1])-ginv(U[,,i])[-(1:2),2])
        
        delta[i] <- 0
        for(j in 1:n){
          delta[i] <- delta[i] + 
            f.z.y[j,i]*
            sum(
              ginv(V[,,i])%*%
                (t(y[,,j]) - matrix(nu[,i], p, r) - matrix(mu[,i], p, r, byrow = TRUE) - beta[,i]%*%t(alfa[,i]))%*%
                ginv(U[,,i])
            )
        }
        delta[i] <- delta[i]/sum(f.z.y[,i])/sum(ginv(U[,,i]))/sum(ginv(V[,,i]))
        
        M[,,i] <- matrix(mu[,i], r, p) + matrix(nu[,i], r, p, byrow = TRUE) + alfa[,i]%*%t(beta[,i]) + delta[i]
      }
      if(model.name[3]=="PI"){
        
        nu[,i] <- 0
        
        for(j in 1:n){
          nu[-1,i] <- nu[-1,i]+
            f.z.y[j,i]* 
            ginv(V[,,i])[-1,]%*%
            (t(y[,,j]) - matrix(mu[,i], p, r, byrow = TRUE) - beta[,i]%*%t(alfa[,i]) - delta[i]
            )%*%
            ginv(U[,,i])%*%
            (rep(1,r))
        }
        
        nu[-1,i] <- ginv(ginv(V[,,i])[-1,-1])%*%nu[-1,i]/sum(f.z.y[,i])/sum(ginv(U[,,i]))
        
        a[,i] <- 0
        
        for(j in 1:n){
          a[,i] <- a[,i]+
            f.z.y[j,i]* 
            t(Q)%*%ginv(U[,,i])%*%
            (y[,,j] - matrix(nu[,i], r, p, byrow = TRUE) - alfa[,i]%*%t(beta[,i]) - delta[i])%*%
            ginv(V[,,i])%*%
            (rep(1,p))
        }
        a[,i] <- ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%(a[,i] -
                                                    ((t(Q[1,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%a[,i])[1,1])/
                                                    ((t(Q[1,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%Q[1,])[1,1])*
                                                    Q[1,])/
          sum(f.z.y[,i])/sum(ginv(V[,,i]))
        mu[-1,i] <- Q[-1,]%*%a[,i] 
        
        beta[,i] <- c(0,1,rep(0,p-2))
        for(j in 1:n){
          beta[-(1:2),i] <- beta[-(1:2),i] + 
            f.z.y[j,i]*ginv(V[,,i])[-(1:2),]%*%
            (t(y[,,j]) - matrix(mu[,i], p, r, byrow = TRUE) - matrix(nu[,i], p, r) - delta[i]
            )%*%
            ginv(U[,,i])%*%
            alfa[,i] 
        }
        beta[-(1:2),i] <-  ginv(ginv(V[,,i])[-(1:2),-(1:2)])%*%
          (beta[-(1:2),i]/sum(f.z.y[,i])/((t(alfa[,i])%*%ginv(U[,,i])%*%alfa[,i])[1,1])-ginv(V[,,i])[-(1:2),2])
        
        b[,i] <- 0
        for(j in 1:n){
          b[,i] <- b[,i] + 
            f.z.y[j,i]*t(Q)%*%ginv(U[,,i])%*%
            (y[,,j] - matrix(mu[,i], r, p) - matrix(nu[,i], r, p, byrow = TRUE) - delta[i]
            )%*%
            ginv(V[,,i])%*%
            beta[,i]
        }
        
        lambda2 <- (sum(f.z.y[,i])*((t(beta[,i])%*%ginv(V[,,i])%*%beta[,i])[1,1]) +
                      ((t(Q[1,] * ((t(Q[2,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%Q[1,])[1,1])/((t(Q[1,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%Q[1,])[1,1]) - Q[2,])%*%
                          ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%b[,i]))[1,1])/
          (((t(Q[2,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%Q[2,])[1,1]) -
             ((t(Q[1,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%Q[2,])[1,1]) *
             ((t(Q[2,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%Q[1,])[1,1]) / 
             ((t(Q[1,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%Q[1,])[1,1]))
        
        lambda1 <- (-((t(Q[1,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%b[,i])[1,1]) -
                      lambda2 * ((t(Q[1,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%Q[2,])[1,1])) /
          ((t(Q[1,])%*%ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%Q[1,])[1,1]) 
        
        b[,i] <-  ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%
          (b[,i] + lambda1 * Q[1,] + lambda2 * Q[2,])/
          sum(f.z.y[,i])/((t(beta[,i])%*%ginv(V[,,i])%*%beta[,i])[1,1])
        
        if(hh == 1) alfa[1:2,i] = 0:1
        alfa[-(1:2),i] <- Q[-(1:2),]%*%b[,i]
        
        delta[i] <- 0
        for(j in 1:n){
          delta[i] <- delta[i] + 
            f.z.y[j,i]*
            sum(
              ginv(V[,,i])%*%
                (t(y[,,j]) - matrix(nu[,i], p, r) - matrix(mu[,i], p, r, byrow = TRUE) - beta[,i]%*%t(alfa[,i]))%*%
                ginv(U[,,i])
            )
        }
        delta[i] <- delta[i]/sum(f.z.y[,i])/sum(ginv(U[,,i]))/sum(ginv(V[,,i]))
        
        M[,,i] <- matrix(mu[,i], r, p) + matrix(nu[,i], r, p, byrow = TRUE) + alfa[,i]%*%t(beta[,i]) + delta[i]
      }
      if(model.name[3]=="PIO"){
        
        nu[,i] <- colSums(matrix(f.z.y[,i],n,p)*apply(y - array(Q%*%b[,i]%*%t(beta[,i]), dim = c(r,p,n)), c(3,2), mean))/sum(f.z.y[,i])
        nu[,i] <- nu[,i] - nu[1,i]
        
        b[,i] <- 0
        for(j in 1:n){
          b[,i] <- b[,i] + 
            f.z.y[j,i]*
            t(Q) %*% ginv(U[,,i]) %*% (y[,,j] - matrix(nu[,i], r, p, byrow = TRUE))  %*% ginv(V[,,i]) %*% beta[,i]
        }
        b[,i] <- ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%b[,i]/((t(beta[,i])%*%ginv(V[,,i])%*%beta[,i])[1,1])/(sum(f.z.y[,i]))
       
        beta[-1,i] <- 0
        for(j in 1:n){
          beta[-1,i] <- beta[-1,i] + 
            f.z.y[j,i]*
            ginv(V[,,i])[-1,] %*% t(y[,,j] - matrix(nu[,i], r, p, byrow = TRUE)) %*% ginv(U[,,i]) %*% Q %*% b[,i]
        }

        beta[-1,i] <-  ginv(ginv(V[,,i])[-1,-1])%*%
          (beta[-1,i]/sum(f.z.y[,i])/((t(Q%*%b[,i])%*%ginv(U[,,i])%*%Q%*%b[,i])[1,1])-ginv(V[,,i])[-1,1])
        
        M[,,i] <- matrix(nu[,i], r, p, byrow = TRUE) + Q%*%b[,i]%*%t(beta[,i])
      }
      if(model.name[3]=="P"){
        P[,,i] <- ginv(t(Q)%*%ginv(U[,,i])%*%Q)%*%t(Q)%*%ginv(U[,,i])%*%
          apply(array(rep(f.z.y[,i], each = r*p), dim = c(r,p,n)) * y, c(1,2), sum)/
          sum(f.z.y[,i])
        M[,,i] <- Q%*%P[,,i]
      }
      tr=sum(diag(V[,,i]))
      V[,,i]=V[,,i]/tr 
      
      
      if (model.name[1]=="GAR") {
        S<-matrix(0,r,r)
        for (j in 1:n) S=S+f.z.y[j,i]*((y[,,j]-M[,,i]))%*%ginv(V[,,i])%*%t(y[,,j]-M[,,i])
        S<-S/(sum(f.z.y[,i]))
        for (h in 2:r) if (h>m) {T[h,(h-m):(h-1),i]<- -ginv(t(S[(h-m):(h-1),(h-m):(h-1)]))%*%S[h,(h-m):(h-1)]} else {T[h,1:(h-1),i]<- -ginv(t(S[1:(h-1),1:(h-1)]))%*%S[h,1:(h-1)]}
        D[,,i]<-diag(diag(T[,,i]%*%S%*%t(T[,,i]))/p)
        U[,,i]=t(T[,,i])%*%ginv(D[,,i])%*%T[,,i]
        U[,,i]=ginv(U[,,i])
      }
      
      if (model.name[1]=="GARI") {
        S<-matrix(0,r,r)
        for (j in 1:n) S=S+f.z.y[j,i]*((y[,,j]-M[,,i]))%*%ginv(V[,,i])%*%t(y[,,j]-M[,,i])
        S<-S/(sum(f.z.y[,i]))
        for (h in 2:r) if (h>m) {T[h,(h-m):(h-1),i]<- -ginv(t(S[(h-m):(h-1),(h-m):(h-1)]))%*%S[h,(h-m):(h-1)]} else {T[h,1:(h-1),i]<- -ginv(t(S[1:(h-1),1:(h-1)]))%*%S[h,1:(h-1)]}
        D[,,i]<-diag(rep(sum(diag(T[,,i]%*%S%*%t(T[,,i])))/(r*p)),r)
        U[,,i]=t(T[,,i])%*%ginv(D[,,i])%*%T[,,i]
        U[,,i]=ginv(U[,,i])
      }
      
      if (model.name[1]=="VVV") { B=0
      for (j in 1:n) B=B+f.z.y[j,i]*((y[,,j]-M[,,i]))%*%ginv(V[,,i])%*%t(y[,,j]-M[,,i])
      U[,,i]<-B/(p*sum(f.z.y[,i]))}
      
      if (model.name[1]=="EEV") { B=0
      for (j in 1:n) B=B+f.z.y[j,i]*((y[,,j]-M[,,i]))%*%ginv(V[,,i])%*%t(y[,,j]-M[,,i])
      DD[,,i]=eigen(B)$vectors
      omega1[,i]=eigen(B)$values}
      
      if (model.name[1]=="VVI") {B=0
      for (j in 1:n) B=B+f.z.y[j,i]*((y[,,j]-M[,,i]))%*%ginv(V[,,i])%*%t(y[,,j]-M[,,i])
      ni=sum(f.z.y[,i])
      A[,,i]=diag(diag(B))/(prod(diag(B))^(1/r))
      lambda[i]=(prod(diag(B))^(1/r))/(p*ni)
      U[,,i]=lambda[i]*A[,,i]}
      
      if (model.name[1]=="VII") {B=0
      for (j in 1:n) B=B+f.z.y[j,i]*((y[,,j]-M[,,i]))%*%ginv(V[,,i])%*%t(y[,,j]-M[,,i])
      ni=sum(f.z.y[,i])
      A[,,i]=diag(r)
      lambda[i]=(sum(diag(B)))/(r*p*ni)
      U[,,i]=lambda[i]*A[,,i]}
      
      
      if (model.name[1]=="EEI" | model.name[1]=="EII" | model.name[1]=="EEE" | model.name[1]=="EGAR" | model.name[1]=="EGARI") { if (i==1) B=0
      for (j in 1:n) B=B+f.z.y[j,i]*((y[,,j]-M[,,i]))%*%ginv(V[,,i])%*%t(y[,,j]-M[,,i])
      }
      
      
      if (model.name[2]=="VVV") { W=0
      for (j in 1:n) W=W+f.z.y[j,i]*(t(y[,,j]-M[,,i]))%*%ginv(U[,,i])%*%(y[,,j]-M[,,i])
      V[,,i]<-W/(r*sum(f.z.y[,i]))}
      
      if (model.name[2]=="EEV") { W=0
      for (j in 1:n) W=W+f.z.y[j,i]*(t(y[,,j]-M[,,i]))%*%ginv(U[,,i])%*%(y[,,j]-M[,,i])
      L[,,i]=eigen(W)$vectors
      omega2[,i]=eigen(W)$values}
      
      if (model.name[2]=="VVI") {W=0
      for (j in 1:n) W=W+f.z.y[j,i]*(t(y[,,j]-M[,,i]))%*%ginv(U[,,i])%*%(y[,,j]-M[,,i])
      ni=sum(f.z.y[,i])
      CC[,,i]=diag(diag(W))/(prod(diag(W))^(1/p))
      csi[i]=(prod(diag(W))^(1/p))/(r*ni)
      V[,,i]=csi[i]*CC[,,i]}
      
      
      if (model.name[2]=="VII") {W=0
      for (j in 1:n) W=W+f.z.y[j,i]*(t(y[,,j]-M[,,i]))%*%ginv(U[,,i])%*%(y[,,j]-M[,,i])
      ni=sum(f.z.y[,i])
      CC[,,i]=diag(p)
      csi[i]=(sum(diag(W)))/(p*r*ni)
      V[,,i]=csi[i]*CC[,,i]}
      
      if (model.name[2]=="EEI" | model.name[2]=="EII" | model.name[2]=="EEE") {if (i==1) W=0
      for (j in 1:n) W=W+f.z.y[j,i]*(t(y[,,j]-M[,,i]))%*%ginv(U[,,i])%*%(y[,,j]-M[,,i])}
      
      w[i]=sum(f.z.y[,i])/n
      
    }
    
    
    for (i in 1:k) {
      
      if (model.name[1]=="EEV") {A[,,i]=diag(rowSums(omega1)/prod(rowSums(omega1))^(1/r))
      lambda[i]= prod(rowSums(omega1))^(1/r)/(n*p)
      U[,,i]=lambda[i]*DD[,,i]%*%A[,,i]%*%t(DD[,,i])}
      
      
      if (model.name[1]=="EEE") U[,,i]=B/(n*p)
      
      if (model.name[1]=="EGAR") {
        S<-B/n
        for (h in 2:r) if (h>m) {T[h,(h-m):(h-1),i]<- -ginv(t(S[(h-m):(h-1),(h-m):(h-1)]))%*%S[h,(h-m):(h-1)]} else {T[h,1:(h-1),i]<- -ginv(t(S[1:(h-1),1:(h-1)]))%*%S[h,1:(h-1)]}
        D[,,i]<-diag(diag(T[,,i]%*%S%*%t(T[,,i]))/p)
        U[,,i]=t(T[,,i])%*%ginv(D[,,i])%*%T[,,i]
        U[,,i]=ginv(U[,,i])
      }
      
      if (model.name[1]=="EGARI") {
        S<-B/n
        for (h in 2:r) if (h>m) {T[h,(h-m):(h-1),i]<- -ginv(t(S[(h-m):(h-1),(h-m):(h-1)]))%*%S[h,(h-m):(h-1)]} else {T[h,1:(h-1),i]<- -ginv(t(S[1:(h-1),1:(h-1)]))%*%S[h,1:(h-1)]}
        D[,,i]<-diag(rep(sum(diag(T[,,i]%*%S%*%t(T[,,i])))/(r*p)),r)
        U[,,i]=t(T[,,i])%*%ginv(D[,,i])%*%T[,,i]
        U[,,i]=ginv(U[,,i])
      }
      
      
      if (model.name[2]=="EEE") V[,,i]=W/(n*r)
      
      if (model.name[1]=="EEI") {A[,,i]=diag(diag(B))/(prod(diag(B))^(1/r))
      lambda[i]=(prod(diag(B))^(1/r))/(p*n)
      U[,,i]=lambda[i]*A[,,i]}
      
      if (model.name[2]=="EEI") {CC[,,i]=diag(diag(W))/(prod(diag(W))^(1/p))
      csi[i]=(prod(diag(W))^(1/p))/(r*n)
      V[,,i]=csi[i]*CC[,,i]}
      
      
      if (model.name[1]=="EII"){A[,,i]=diag(r)
      lambda[i]=(sum(diag(B)))/(r*p*n)
      U[,,i]=lambda[i]*A[,,i]}
      
      if (model.name[2]=="EEV") {CC[,,i]=diag(rowSums(omega2)/prod(rowSums(omega2))^(1/p))
      csi[i]= prod(rowSums(omega2))^(1/p)/(n*r)
      V[,,i]=csi[i]*L[,,i]%*%CC[,,i]%*%t(L[,,i])}
      
      
      if (model.name[2]=="EII"){CC[,,i]=diag(p)
      csi[i]=(sum(diag(W)))/(p*r*n)
      V[,,i]=csi[i]*CC[,,i]}
      
    }
    for (i in 1:k) f.y.z[,i]=dmm(y,M[,,i],U[,,i],V[,,i], log = TRUE)
    f.y.z=ifelse(is.na(f.y.z),mean(f.y.z, na.rm=T),f.y.z)
    
    temp = sum(rowLogSumExps(matrix(log(w),n,k,byrow=TRUE) + f.y.z)) #sum(log(w%*%t(f.y.z)))
    if(temp == - Inf){
      temp = -.Machine$double.xmax
      #cat(paste0("iter ", hh, ": temp == - Inf \n"))
    }
    likelihood<-c(likelihood,temp)
    ratio<-abs((temp-lik)/lik)
    #if (temp < lik) cat(paste0("iter ", hh, ": temp < lik \n"))
    #if ((temp < lik) & (hh > 7)) ratio<-eps
    if(temp > lik) lik<-temp
    if(hh == it) if(likelihood[hh] < likelihood[hh-1])  it <- it+1
  }
  
  w=w.old
  V=V.old
  U=U.old
  M=M.old
  f.z.y=f.z.y.old
  
  lik = sum(rowLogSumExps(matrix(log(w),n,k,byrow=TRUE) + f.y.z)) #sum(log(w%*%t(f.y.z)))
  
  #if(lik == - Inf) lik = -.Machine$double.xmax
  if(!is.null(test)){
    f.y.z.test <- matrix(0,dim(y.test)[3],k)
    for (i in 1:k) f.y.z.test[,i]=dmm(y.test,M[,,i],U[,,i],V[,,i], log = TRUE)
    f.y.z.test=ifelse(is.na(f.y.z.test),mean(f.y.z.test, na.rm=T),f.y.z.test)
    lik.test <- sum(rowLogSumExps(matrix(log(w),dim(y.test)[3],k,byrow=TRUE) + f.y.z.test))
  }else{lik.test <- NULL}
  
  conta.zeri<-0
  for (h in 1:r) if ((h+m+1)<=r) conta.zeri<-conta.zeri+(r-h-m)
  
  if (model.name[1]=="GAR")  h1=k*r+k*(r*(r-1)/2)-k*conta.zeri
  if (model.name[1]=="GARI") h1=k+k*(r*(r-1)/2)-k*conta.zeri
  if (model.name[1]=="EGAR")  h1=r+(r*(r-1)/2)-conta.zeri
  if (model.name[1]=="EGARI") h1=1+(r*(r-1)/2)-conta.zeri
  if (model.name[1]=="III")  h1=0
  if (model.name[2]=="III")  h2=0
  if (model.name[1]=="EEV")  h1=(k*r*(r+1)/2)-((k-1)*r)
  if (model.name[2]=="EEV")  h2=(k*p*(p+1)/2)-((k-1)*p)
  if (model.name[1]=="VVV")  h1=(k*r*(r+1)/2)
  if (model.name[2]=="VVV")  h2=(k*p*(p+1)/2)
  if (model.name[1]=="EEE")  h1=(r*(r+1)/2)
  if (model.name[2]=="EEE")  h2=(p*(p+1)/2)
  if (model.name[1]=="VVI")  h1=(k*r)
  if (model.name[2]=="VVI")  h2=(k*p)
  if (model.name[1]=="EEI")  h1=r
  if (model.name[2]=="EEI")  h2=p
  if (model.name[1]=="VII")  h1=k
  if (model.name[2]=="VII")  h2=k
  if (model.name[1]=="EII")  h1=1
  if (model.name[2]=="EII")  h2=1
  if (model.name[3]=="U")  h3=p*r
  if (model.name[3]=="P")  h3=p*q
  if (model.name[3]=="A" | model.name[3]=="I")  h3=p+r-1
  if (model.name[3]=="GI")  h3=2*p+2*r-5
  if (model.name[3]=="PI")  h3=2*p+2*q-5
  if (model.name[3]=="PIO")  h3=2*p+q-2
  
  #h=k-1+k*(p*r)+h1+h2
  h=k-1+k*h3+h1+h2
  
  aic=-2*lik+2*h
  bic=-2*lik+h*log(n)
  EN=entr(f.z.y)
  icl.bic=-2*lik+2*EN+h*log(n)
  cl=apply(f.z.y,1,which.max)
  output<-list(likelihood=likelihood,lik.test=lik.test,w=w,U=U,nu=nu,mu=mu,beta=beta,alfa=alfa,delta=delta,a=a,b=b,P=P,M=M,V=V,f.z.y=f.z.y,cl=cl,bic=bic,aic=aic,icl.bic=icl.bic,tempo=proc.time()-ptm,D=D,T=T)
  
  invisible(output)
  
}


array.matrix<-function(A,B)
{
  ## A is an array of dimensions k1 x k2 x k3
  ## B is a matrix of dimensions k3 x k4
  ## returns D of dimensions k1 x k2 x k4
  
  C<-A%o%B
  D<-apply(C,c(1,2,5),diag)
  D<-(apply(D,c(2,3,4),sum))
  return(D)
}



rmm=function(r,M,U,V)
{
  n=nrow(M)
  p=ncol(M)
  z=array(rnorm(n*p*r),c(n,p,r))
  y=array(0,c(n,p,r))
  for (i in 1:r) {
    y[,,i]=t(chol(U))%*%z[,,i]%*%(chol(V))+M
  }
  return(y)
}

### NOTA y è di dimensione r x p x numobs per ottenere la varianza U%x%V si fa cosi:
### y=rmm(10000,M,U,V)
### y=aperm(y,c(2,1,3))
### yy=t(apply(y,3,c))
### var(yy)


dmm=function(y,M,U,V,log=FALSE)
{
  #if (det(V)<1.0e-30){ 
  #print(V)
  #  diag(V)=diag(V)+0.1
  #  cat(paste0("det(V)<1.0e-30 \n"))
  #}
  #  if(!isSymmetric(V)){
  #    cat(paste0("!isSymmetric(V) \n"))
  #    V[lower.tri(V)] <- t(V)[lower.tri(V)]
  #  }
  
  #if (det(U)<1.0e-30){
  #print(U)
  #  diag(U)=diag(U)+0.1
  #  cat(paste0("det(U)<1.0e-30 \n"))
  #}
  
  # if(!isSymmetric(U)){
  #    cat(paste0("!isSymmetric(U) \n"))
  #    U[lower.tri(U)] <- t(U)[lower.tri(U)]
  #  }
  #if(!isSymmetric(U)){ 
  #    print(U)
  # }
  
  r=dim(y)[1]
  p=dim(y)[2]
  n=dim(y)[3]
  #if(!isSymmetric(U)) 
  #print(U)
  #if(!isSymmetric(V)) 
  #print(V)
  sigma=V%x%U
  if (det(sigma)<1.0e-30){
    #cat(paste0(isSymmetric(sigma), " det(V%x%U)<1.0e-30 \n"))
    diag(sigma)=diag(sigma)+0.1
  }
  if(!isSymmetric(sigma)){
    sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
    #if(!isSymmetric(sigma)) {print(sigma)
    #for(i in 1:(p*r - 1)){
    #  for(j in (i+1):(p*r)){
    #   if(sigma[i,j] != sigma[j,i]) cat(paste(i, j, sigma[i,j], sigma[j,i], "\n"))
    #  }
    #}
    #}
  }
  y=t(apply(y,3,c))
  mu=c(M)
  if (log) f.y=dmvnorm(y,mu,sigma,log=TRUE) else f.y=dmvnorm(y,mu,sigma)
  #if (log) f.y=dmatrixnorm(y,M,U = U, V = V, log=TRUE) else f.y=dmatrixnorm(y,M,U = U, V = V)
  return(f.y)
}

misc=function(classification, truth)
{
  q <- function(map, len, x) {
    x <- as.character(x)
    map <- lapply(map, as.character)
    y <- sapply(map, function(x) x[1])
    best <- y != x
    if (all(len) == 1)
      return(best)
    errmin <- sum(as.numeric(best))
    z <- sapply(map, function(x) x[length(x)])
    mask <- len != 1
    counter <- rep(0, length(len))
    k <- sum(as.numeric(mask))
    j <- 0
    while (y != z) {
      i <- k - j
      m <- mask[i]
      counter[m] <- (counter[m]%%len[m]) + 1
      y[x == names(map)[m]] <- map[[m]][counter[m]]
      temp <- y != x
      err <- sum(as.numeric(temp))
      if (err < errmin) {
        errmin <- err
        best <- temp
      }
      j <- (j + 1)%%k
    }
    best
  }
  if (any(isNA <- is.na(classification))) {
    classification <- as.character(classification)
    nachar <- paste(unique(classification[!isNA]), collapse = "")
    classification[isNA] <- nachar
  }
  MAP <- mapClass(classification, truth)
  len <- sapply(MAP[[1]], length)
  if (all(len) == 1) {
    CtoT <- unlist(MAP[[1]])
    I <- match(as.character(classification), names(CtoT),
               nomatch = 0)
    one <- CtoT[I] != truth
  }
  else {
    one <- q(MAP[[1]], len, truth)
  }
  len <- sapply(MAP[[2]], length)
  if (all(len) == 1) {
    TtoC <- unlist(MAP[[2]])
    I <- match(as.character(truth), names(TtoC), nomatch = 0)
    two <- TtoC[I] != classification
  }
  else {
    two <- q(MAP[[2]], len, classification)
  }
  err <- if (sum(as.numeric(one)) > sum(as.numeric(two)))
    as.vector(one)
  else as.vector(two)
  bad <- seq(along = classification)[err]
  errorRate = length(bad)/length(truth)
  return(errorRate)
}


mapClass=function (a, b)
{
  l <- length(a)
  x <- y <- rep(NA, l)
  if (l != length(b)) {
    warning("unequal lengths")
    return(x)
  }
  aChar <- as.character(a)
  bChar <- as.character(b)
  Tab <- table(a, b)
  Ua <- dimnames(Tab)[[1]]
  Ub <- dimnames(Tab)[[2]]
  aTOb <- rep(list(Ub), length(Ua))
  names(aTOb) <- Ua
  bTOa <- rep(list(Ua), length(Ub))
  names(bTOa) <- Ub
  k <- nrow(Tab)
  Map <- rep(0, k)
  Max <- apply(Tab, 1, max)
  for (i in 1:k) {
    I <- match(Max[i], Tab[i, ], nomatch = 0)
    aTOb[[i]] <- Ub[I]
  }
  if (is.numeric(b))
    aTOb <- lapply(aTOb, as.numeric)
  k <- ncol(Tab)
  Map <- rep(0, k)
  Max <- apply(Tab, 2, max)
  for (j in (1:k)) {
    J <- match(Max[j], Tab[, j])
    bTOa[[j]] <- Ua[J]
  }
  if (is.numeric(a))
    bTOa <- lapply(bTOa, as.numeric)
  list(aTOb = aTOb, bTOa = bTOa)
}


entr<-function(z)
{
  numobs<-nrow(z)
  numg<-ncol(z)
  temp<-0
  z<-ifelse(z==0,z+0.000000000000000000000001,z)
  for (i in 1:numg) for (j in 1:numobs) temp<-temp+(z[j,i]*log(z[j,i]))
  return(-temp)
}
