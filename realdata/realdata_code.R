rm(list=ls(all=TRUE))
library(countyweather)
library(usmap)
library(stringr)

load("countydata0401.rdata")
fipscd = rep(NA, nrow(county))
county$fips= str_pad(county$fips, 5, pad = '0')
fipscd = county$fips


load("marchdaily1.rdata")

###process the data after load daily data
rlnames = unique(dailydata[, 2])
clnames = unique(dailydata[, 3])
dataprc = matrix(NA, length(rlnames), length(clnames))
dataprc = data.frame(dataprc,  row.names = rlnames)
colnames(dataprc) = clnames
datatmax = datatmin = dataprc
for(i in 1:nrow(dailydata)){
    x = as.character(dailydata[i, 2])
    y = as.character(dailydata[i, 3])
    dataprc[x, y] = dailydata[i, 4]
    datatmax[x, y] = dailydata[i, 5]
    datatmin[x, y] = dailydata[i, 6]
}

unqfips = dailydata[!duplicated(dailydata[, 2]), 1]
dataprc = cbind(unqfips, dataprc)
datatmax = cbind(unqfips, datatmax)
datatmin = cbind(unqfips, datatmin)
ix = which(apply(is.na(dataprc[, -1]), 1, mean) >0.6 |apply(is.na(datatmax[, -1]), 1, mean) >0.6 |apply(is.na(datatmin[, -1]), 1, mean) >0.6 )
dataprc = dataprc[-ix, ]
datatmax = datatmax[-ix, ]
datatmin = datatmin[-ix, ]
mapprox= function(x){
    temp = approxfun(1:length(x), x, rule = 2)(1:length(x))
    return(temp)
    }
dataprc[, -1] = t(apply(dataprc[, -1], 1, mapprox))
datatmax[, -1] = t(apply(datatmax[, -1], 1, mapprox))
datatmin[, -1] = t(apply(datatmin[, -1], 1, mapprox))
datadiff = datatmax[, -1] - datatmin[, -1]
dataavg = datatmax[, -1]/2 +  datatmin[, -1]/2

repdesign = cbind(dataprc, datadiff, dataavg)
unqfips = unique(repdesign[, 1])
mdesign = matrix(NA, length(unqfips), ncol(repdesign[, -1]))

omega = array(NA, c(ncol(repdesign[, -1]), ncol(repdesign[, -1]), length(unqfips)))
for(i in 1:length(unqfips)){
    ix = repdesign[, 1] == unqfips[i]
    omega[, , i] = cov(repdesign[ix, -1])/(sum(ix)-1)
    mdesign[i, ]= apply(repdesign[ix, -1, drop = FALSE], 2, mean)
    
}
momega = apply(omega, c(1, 2), mean,trim = 0.2,  na.rm = T)
momega = apply(omega, c(1, 2), median, na.rm = T)

selected_Index<-c()
for (i in 1:length(unqfips)){
  Index_update<-which(county$fips==unqfips[i])
  selected_Index<-rbind(selected_Index,Index_update)
}

resp = county[selected_Index, ]



###### 1-31  precipitation
###### 32-62 daily temperature difference
#####  63-93 daily temperature average



#library(ncvreg)
library(glmnet)
library(mvtnorm)
#library(parallel)
library(MASS)
source("simplex.R")
source("functionserroroff.R")

##infected cases
mytol = 1e-5
maxitr = 300
gamma = 1



lamlist = seq(0.01, 0.05, length = 5)/10
lamlist1 = seq(0.1, 0.5, length = 5)

Y=resp$cases
n=length(Y)
X=cbind(1,  mdesign)
p=dim(X)[2]
offset=log(resp$POPESTIMATE2019)


tempm = apply(omega, c(1, 2), median, na.rm = T)
temp = eigen(tempm)
temp$values[temp$values<0] = 0
tempm = temp$vectors %*% diag(temp$values) %*% solve(temp$vectors)
momega = cbind(0, tempm)
momega = rbind(0, momega)

rep_t=120
set.seed(1000)
TTres = array(NA, dim=c(rep_t, 94, 2))
TTres0 = array(NA, dim=c(rep_t, 94, 2))

for (i in 1:rep_t){
  print(i)
  ##add error
  U<-rmvnorm(n=n,  mean=rep(0,94), sigma=momega*0.9)
  W<-X+U
  
  
  maxoff = offset#max(offset)
  for(j in 2:94){
    Index_N<- j
    fit<-glmnet(W, Y/exp(maxoff),  alpha = 1,offset = offset - maxoff,  lambda =min(lamlist), intercept = F, family = 'poisson')#glm(resp$cases~mdesign + offset(log(resp$POPESTIMATE2019)) -1,   family = 'quasipoisson')
    beta00<-as.matrix(coef(fit))[-1]
    beta00 <-euclidean_proj_l1ball(beta00,s=10*sqrt(2))
    beta00=shrink(beta00,10)	
    fit<-glmnet(W[, -Index_N], Y/exp(maxoff), offset = offset - maxoff,  alpha = 1, lambda = min(lamlist1), intercept = F, family = 'poisson')#glm(resp$cases~mdesign + offset(log(resp$POPESTIMATE2019)) -1,   family = 'quasipoisson')
    beta01<-as.matrix(coef(fit))[-1]
    beta01 <-euclidean_proj_l1ball(beta01,s=10*sqrt(2))
    beta01=shrink(beta01,10)
  
    
    C<- matrix(0,1,p)
    C[,Index_N]=1
    b <-0
    N <-  c(rep(FALSE, (Index_N-1)),TRUE, rep(FALSE, (p-Index_N)))
    R <- 100 * norm(beta01,type='2')
    N1 <- c(rep(FALSE, (Index_N-1)),FALSE, rep(FALSE, (p-Index_N)))
    mmomega = momega[, -Index_N]
    mmomega = mmomega[-Index_N, ]
    
    beta.uncon<-try(cv.SCAD_ADMM_uncon(W=W, Y=Y/exp(maxoff), offset=offset- maxoff,N=N,omega=momega, beta0=beta00, R=R,err=1e-4, tune="BIC", lambdalist =lamlist1, rho = 0.1))
    beta.con<-try(cv.SCAD_ADMM_uncon(W=W[, -Index_N], Y=Y/exp(maxoff), offset=offset - maxoff, N=N[-Index_N], omega=mmomega, beta0=beta01, R=R,err=1e-4, tune="BIC", lambdalist= lamlist1, rho = 0.1))
    if (class(beta.uncon)=="try-error"|class(beta.con)=="try-error"){next}
    else{
      if(Index_N<p){
      beta.con = c(beta.con[1:(Index_N-1)], 0, beta.con[Index_N:93])
    }else{beta.con = c(beta.con[1:(Index_N-1)], 0)}
    #W=X
    index.uncon<- which(beta.uncon!=0)
    index.con<- which(beta.con!=0)
    
    ##Wald
    if(sum(index.uncon)>0){
      Qwald<-Qest(W, offset - maxoff,momega, beta.uncon)
      Sigwald<-Sigest(W, Y/exp(maxoff), offset - maxoff, momega, beta.uncon)
      invPhi<-try(solve(Phiest(Sigwald,Qwald,beta.uncon,C,Index_N)))
      if(inherits(invPhi, "try-error")){
        invPhi <- ginv(Phiest(Sigwald,Qwald,beta.uncon,C,Index_N))
      }
      TW<-n*t(C%*%beta.uncon-b)%*%invPhi%*%(C%*%beta.uncon-b)
    }else{TW<-0}
    
    ##Score
    DL <- DL1(W, Y/exp(maxoff), offset - maxoff,momega, beta.con)
    index <- union(index.con,Index_N)
    
    if(sum(index.con)>0){
      Sigscore<-Sigest(W, Y/exp(maxoff), offset - maxoff,momega, beta.con)
      Qscore<-Qest(W,offset-maxoff, momega, beta.con)
      invPhi<-try(solve(Phiest(Sigscore,Qscore,beta.con,C,Index_N)))
      L = C[, index] %*% solve(Qscore[index, index])
      midd = t(L) %*% invPhi %*% L
      TS<-n*DL[index]%*%midd%*%DL[index]
    }else{TS<-0}
    print(c(TS, TW))
    TTres[i,Index_N, ]= c(TS, TW)}
  }
  
  ###Without error
  
  for(j in 2:94){
    Index_N<- j
    fit<-glmnet(W, Y/exp(maxoff),  offset = offset - maxoff, alpha = 1, lambda =min(lamlist1), intercept = F, family = 'poisson')#glm(resp$cases~mdesign + offset(log(resp$POPESTIMATE2019)) -1,   family = 'quasipoisson')
    beta00<-as.matrix(coef(fit))[-1]
    beta00 <-euclidean_proj_l1ball(beta00,s=10*sqrt(2))
    beta00=shrink(beta00,10)	
    fit<-glmnet(W[, -Index_N], Y/exp(maxoff), offset = offset - maxoff,  alpha = 1, lambda = min(lamlist), intercept = F, family = 'poisson')#glm(resp$cases~mdesign + offset(log(resp$POPESTIMATE2019)) -1,   family = 'quasipoisson')
    beta01<-as.matrix(coef(fit))[-1]
    beta01 <-euclidean_proj_l1ball(beta01,s=10*sqrt(2))
    beta01=shrink(beta01,10)	
    momega0=matrix(0,p,p)
    mmomega0 = momega0[, -Index_N]
    mmomega0 = mmomega0[-Index_N, ]
    
    C<- matrix(0,1,p)
    C[,Index_N]=1
    b <-0
    N <-  c(rep(FALSE, (Index_N-1)),TRUE, rep(FALSE, (p-Index_N)))
    R <- 100*norm(beta00,type='2')
    
    
    
    
    beta.uncon0<-try(cv.SCAD_ADMM_uncon(W=W, Y=Y/exp(maxoff), offset=offset - maxoff,N=N,omega=momega0, beta0=beta00, R=R,err=1e-4, tune="BIC", lambdalist = lamlist1, rho = 0.1))
    beta.con0<-try(cv.SCAD_ADMM_uncon(W=W[, -Index_N], Y=Y/exp(maxoff), offset=offset -maxoff,N=N[-Index_N],omega=mmomega0, beta0=beta01, R=R,err=1e-4, tune="BIC", lambdalist = lamlist1, rho = 0.1))
    if (class(beta.uncon0)=="try-error"|class(beta.con0)=="try-error"){next}
    else{
    if(Index_N <94){
      beta.con0 = c(beta.con0[1:(Index_N-1)], 0, beta.con0[Index_N:93])
    }else{
      beta.con0 = c(beta.con0[1:(Index_N-1)], 0)
    }
    #W=X
    index.uncon<- which(beta.uncon0!=0)
    index.con<- which(beta.con0!=0)
    ##Wald
    if(sum(index.uncon)>0){
      Qwald<-Qest(W, offset - maxoff,momega0, beta.uncon0)
      Sigwald<-Sigest(W, Y/exp(maxoff), offset - maxoff, momega0, beta.uncon0)
      invPhi<-try(solve(Phiest(Sigwald,Qwald,beta.uncon0,C,Index_N)))
      if(inherits(invPhi, "try-error")){
        invPhi <- ginv(Phiest(Sigwald,Qwald,beta.uncon0,C,Index_N))
      }
      TW<-n*t(C%*%beta.uncon0-b)%*%invPhi%*%(C%*%beta.uncon0-b)
    }else{TW<-0}
    
    
    
    DL <- DL1(W, Y/exp(maxoff), offset - maxoff,momega0, beta.con0)
    index <- union(index.con,Index_N)
    
    if(sum(index.con)>0){
      Sigscore<-Sigest(W, Y/exp(maxoff), offset-maxoff,momega0, beta.con0)
      Qscore<-Qest(W,offset-maxoff, momega0, beta.con0)
      invPhi<-try(ginv(Phiest(Sigscore,Qscore,beta.con0,C,Index_N)))
      L = C[, index] %*% ginv(Qscore[index, index])
      midd = t(L) %*% invPhi %*% L
      TS<-n*DL[index]%*%midd%*%DL[index]
    }else{TS<-0}
    print(c(TS, TW))
    TTres0[i,Index_N, ]= c(TS, TW)}
  } 
}
