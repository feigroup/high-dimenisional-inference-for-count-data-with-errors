rm(list=ls(all=TRUE))
library(mvtnorm)
#library(ncvreg)
library(MASS)
source("simplex.R")
source("functionserror.R")
set.seed(2018)
p<-50
n<-300
cor<-0.5
rep<-500
omega<-matrix(0,p,p)
diag(omega) = 0.05
beta_con<-matrix(0,rep,p)
beta_uncon<-matrix(0,rep,p)
TW_all<-matrix(0,rep,1)
TS_all<-matrix(0,rep,1)

for (ii in 1:rep){
  print(ii)
  U<-rmvnorm(n=n,  mean=rep(0,p), sigma=omega)
  X<-rmvnorm(n=n,  mean=rep(0,p), sigma=diag(0.5,p,p))
  W<-X + U
  beta_true<-matrix(rep(0,p),p,1)
  beta_true[1,]<- 0.75
  beta_true[2,]<- -0.75

  posmean<-X%*%beta_true
  Y<-rpois(n, exp(posmean))
  
  #Test beta_1+beta_2=0
  C <- t(as.matrix(c(1,1, rep(0,p-2))))
  b <- 0
  N <- c(TRUE, TRUE, rep(FALSE, p-2))
  R <- 1.5*norm(beta_true,type='2')
  Index_N<- c(1,2)
  
  
  beta.uncon<-cv.SCAD_ADMM_uncon(W=W, Y=Y, N=N,omega=omega, beta0=rep(0, dim(W)[2]), R=R,err=1e-4, tune="BIC")
  beta_uncon[ii,]<- beta.uncon
  index.uncon<- which(beta.uncon!=0)
  beta.con=cv.SCAD_ADMM_con(W=W, Y=Y, N=N, omega=omega, C=C, b=b, beta0=rep(0, dim(W)[2]), R=R,err=1e-4, tune="BIC")
  beta_con[ii,]<- beta.con
  index.con<- which(beta.con!=0)
  ###construct test 
  ##Wald 
  if(sum(index.uncon)>0){
    Qwald<-Qest(W, omega, beta.uncon)
    Sigwald<-Sigest(W, Y, omega, beta.uncon)
    invPhi<-try(solve(Phiest(Sigwald,Qwald,beta.uncon,C,Index_N)))
    if(inherits(invPhi, "try-error")){
      invPhi <- INV(Phiest(Sigwald,Qwald,beta.uncon,C,Index_N))
    }
    TW<-n*t(C%*%beta.uncon-b)%*%invPhi%*%(C%*%beta.uncon-b)
  }else{TW<-0}
  TW_all[ii,]<-TW
  #print(TW)
  ##Score
  DL <- DL1(W, Y, omega, beta.con)
  index <- union(index.con,Index_N)
  
  if(sum(index.con)>0){
      Sigscore<-Sigest(W, Y, omega, beta.con)
      Qscore<-Qest(W, omega, beta.con)
      invPhi<-try(solve(Phiest(Sigscore,Qscore,beta.con,C,Index_N)))
      L = C[, index] %*% solve(Qscore[index, index])
      midd = t(L) %*% invPhi %*% L
    TS<-n*DL[index]%*%midd%*%DL[index]
  }else{TS<-0}
  
  TS_all[ii,]<-TS
  #print(TS)
}
cvalue=qchisq(0.95,df=1)
mean(TW_all>cvalue)
mean(TS_all>cvalue)
