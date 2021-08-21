source("simplex.R")
####################################################################################
## Estimate parameters using ADMM for high-dimensional poisson regression 
## with the linear constraints: C beta=b
####################################################################################
## Objective function:
## f(Y,X,beta)+p_a(lambda,theta_N^c)+rho v_1^T (beta-theta)+rho v_2^T (Cbeta-b)+0.5rho ||beta-theta||_2^2
## subject to C beta=b
## f: the loglikelihood
## p: The SCAD penalty function 
## a: parameter in the SCAD penalty, default=3.7
## Y: the response vector
## X: the design matrix
## lambda: the regularization parameter
## rho: regularization parameter in ADMM, we set rho=1
## v_1,v_2: the Langrange multiplier
## #N={1,2}
#####################################################################################
# ADMM step: update beta

## Inversion of a matrix
INV <- function(w) {as.matrix(t(svd(w)$v %*% (t(svd(w)$u)/svd(w)$d)))}


## Derivative of likelihood
# DL_function<-function(Y,W,omega,beta){
#   first_term=crossprod(W/n, Y)
#   exp_term1=exp(W%*%beta-as.numeric(t(beta)%*%omega%*%beta/2))
#   exp_term2=W-t(matrix(rep(t(beta)%*%omega,n),p,n))
#   exp_term=crossprod(exp_term2/n, exp_term1)
#   ##p-by-1
#   DL1_likelihood=-(first_term-exp_term)
#   
#   exp_term_DL21=t(exp_term2)%*%diag(c(exp_term1))%*%exp_term2
#   exp_term_DL22=sum(exp_term1)*omega
#   ##p-by-p
#   DL2_likelihood=(exp_term_DL21-exp_term_DL22)/n
#   
#   DL_likelihood=list("","")
#   DL_likelihood[[1]]=DL1_likelihood
#   DL_likelihood[[2]]=DL2_likelihood
#   return(DL_likelihood)
# }


##update_beta
one_step_update_beta_con <- function(W, Y, offset, omega, beta, theta, C, b, v1, v2, rho=1){
  n <- length(Y)
  p <- dim(W)[2] 
  
  first_term=crossprod(W/n, Y)
  exp_term1=exp(W%*%beta+offset-as.numeric(t(beta)%*%omega%*%beta/2))
  exp_term2=W-t(matrix(rep(t(beta)%*%omega,n),p,n))
  exp_term=crossprod(exp_term2/n, exp_term1)
  ##p-by-1
  DL1_likelihood=-(first_term-exp_term)
  
  exp_term_DL21=t(exp_term2)%*%diag(c(exp_term1))%*%exp_term2
  exp_term_DL22=sum(exp_term1)*omega
  ##p-by-p
  DL2_likelihood=(exp_term_DL21-exp_term_DL22)/n
  
  beta_new <- (rho*theta-v1-DL1_likelihood+as.vector(t(C)%*%(rho * b-v2)) + DL2_likelihood%*%beta)	   
  if (p<n){
    #beta_new <- solve( crossprod(X,as.vector(lambda)*X)/n+rho*diag(p), beta_new/2)
    B0 <- DL2_likelihood+rho*diag(p)+rho*crossprod(C,C)
    beta_new1 = ginv(t(B0) %*% B0, tol = mytol) %*% t(B0) %*% (beta_new)
    #beta_new1 <- ginv(B0, tol = 1e-5)%*%beta_new
      #try(solve(B0, beta_new))
    #if(inherits(beta_new1, "try-error"))
    #  beta_new1 <- ginv(DL2_likelihood+rho*diag(p)+rho*crossprod(C,C))%*%beta_new
    beta_new <- beta_new1
  }
  else{
    # using woodbury formula when p is large
    A<-rho*diag(p)-omega*mean(exp_term1)+rho*crossprod(C,C)
    U<-t(W-t(matrix(rep(t(beta)%*%omega,n),p,n)))
    C0<-diag(n)/n
    V<-as.vector(exp_term1)*(W-t(matrix(rep(t(beta)%*%omega,n),p,n)))
    invA<-ginv(A,tol=0.0001)
    #try(solve(A))
    #if(inherits(invA, "try-error")){
    #  invA <- ginv(A)
    #}
    invmain<-ginv(n*diag(n)+V%*%invA%*%U,tol=0.0001)
    #try(solve(n*diag(n)+V%*%invA%*%U))
    #if(inherits(invmain, "try-error")){
    #  invmain <- ginv(n*diag(n)+V%*%invA%*%U)
    #}
    beta_new<-(invA-invA%*%U%*%invmain%*%V%*%invA)%*%beta_new
  }
  
  return (beta_new)
}

update_beta_con <- function(W, Y, offset,beta, omega,theta, C, b, v1, v2, rho=1){
  
  diff <- 1
  count <- 0
  
  while(diff>1e-4&&count<=200){
    beta_new <- one_step_update_beta_con(W=W, Y=Y, offset=offset,omega=omega, beta=beta, theta=theta,  C=C, b=b, v1=v1, v2=v2, rho=1)
    diff <- sqrt(sum((beta_new-beta)^2))
    #    print(diff)
    if (is.na(diff)||is.nan(diff)||is.infinite(diff))
      return(rep(0, length(beta)))
    beta <- beta_new
    count <- count+1
  }
  
  return (beta)
}

# ADMM step: update theta
# SCAD thresholding
update_theta_con <- function(beta, v1, N, lambda, rho=1, a=3.7){
  theta <- v1/rho+beta * rho
  #N <- c(TRUE, TRUE, rep(FALSE, length(theta)-2))
  # when absolute value smaller than lambda
  theta[(!N) & abs(theta)<=lambda] <- 0
  # when larger than lambda but smaller than 2lambda
  index <- (!N) & abs(theta)>lambda & abs(theta)<=2*lambda
  theta[index] <- sign(theta[index])*(abs(theta[index])-lambda)
  # when larger than 2lambda but smaller than a lambda
  index <- (!N) & abs(theta)>2*lambda & abs(theta)<=a*lambda
  theta[index] <- ((a-1)*theta[index]-sign(theta[index])*a*lambda)/(a-2) 
  
  return (theta)
}

# calculate the primal and dual residuals
res_con <- function(beta_new, theta_new, beta_old, theta_old, C, b, rho=1){
  # the primal residual
  r <- c(beta_new-theta_new)
  # the dual residual 
  s <- rho*(theta_new-theta_old)
  
  return (list(r=r, s=s))
}

# ADMM algorithm
SCAD_ADMM_con <- function(W, Y, offset, N, omega, C, b, beta0, lambda, R, rho=1, a=3.7, err=0.5e-3){
  
  # initial value of parameters
  beta_old <- beta0
  theta_old <- beta0
  v1 <- 0
  v2 <- 0
  iter <- 0
  
  eps <- list(r=1, s=1)
  while ((norm(as.matrix(eps$r), "F")>=err ||  (norm(as.matrix(eps$s), "F")>=err))&&iter<=maxitr)#
  {
    # ADMM update
    # update v
    v1 <- v1+gamma *  (beta_old-theta_old)
    v2 <- v2+gamma *  (C%*%beta_old-b)
    
    # update beta
    beta_new <- update_beta_con(W=W, Y=Y, offset=offset, omega=omega, beta=beta_old, theta=theta_old, C=C, b=b, v1=v1, v2=v2, rho=rho)
    
    # update theta
    theta_new <- update_theta_con(beta=beta_new, v1=v1, N=N, lambda=lambda, rho=rho, a=a)
    
    #project
    theta_new<-euclidean_proj_l1ball(theta_new,s=R*sqrt(2))
    theta_new=shrink(theta_new,R)
    
    # calculate the residuals
    eps <- res_con(beta_new, theta_new, beta_old, theta_old, C=C, b=b, rho=rho)
    #print(norm(as.matrix(eps$r), "F"))
    #print(norm(as.matrix(eps$s), "F"))
   # print(theta_new)
    # beta_old and theta_old
    beta_old <- beta_new
    theta_old <- theta_new
    
    iter <- iter+1
    #print(iter)
  }
  return(theta_new)
}
##############################################################################
## K-folded cross validation, K=5 by default
##############################################################################
cv.SCAD_ADMM_con <- function(W, Y, offset, N, omega, C, b, beta0, R, K=5, rho=1, a=3.7, err=0.5e-3, tune="BIC", lambdalist){
  
  # potential tuning parameters
  #lambdalist <- exp(seq(-2.5, 0.5, 0.075))
  lambdan <- length(lambdalist)
  
  # data splitting
  n <- length(Y)
  p <- dim(W)[2]
  
  beta.all <- matrix(0, lambdan, p)
  
  if (tune=="cv"){
    folds <- split(sample(n, n, replace=FALSE), as.factor(1:K))
    
    # calculate MSE for each folds
    MSE <- rep(0, lambdan)
    for (k in 1:K){
      # testing dataset
      W0 <- W[folds[[k]], ]
      Y0 <- Y[folds[[k]]]
      offset0<-offset[folds[[k]]]
      
      # training dataset
      W1 <- W[-folds[[k]], ]
      Y1 <- Y[-folds[[k]]]
      offset1<-offset[-folds[[k]]]
      
      # est the training dataset
      for (j in 1:lambdan){
        beta0 <- SCAD_ADMM_con(W=W1, Y=Y1, offset=offset1, N=N, omega=omega, C=C, b=b, beta0=beta0, lambda=lambdalist[j], R=R,
                              rho=rho, a=a, err=err)
        MSE[j] <- MSE[j]+sum(exp(W%*%beta0+offset-as.numeric(t(beta0)%*%omega%*%beta0/2))-Y*(W%*%beta0+offset))
      }
    }
    
    # take the minimum lambda
    lambda <- lambdalist[which.min(MSE)]
    beta <- SCAD_ADMM_con(W=W, Y=Y, offset=offset, N=N, omega=omega, C=C, b=b, beta0=beta0, lambda=lambda, R=R, rho=rho, a=a, 
                         err=err)
  }
  else{
    BIC <- rep(0, lambdan)
    for (j in 1:lambdan){
      beta.all[j,] <- SCAD_ADMM_con(W=W, Y=Y, offset=offset,N=N, omega=omega, C=C, b=b, beta0=beta0, lambda=lambdalist[j], R=R,
                                   rho=rho, a=a, err=err)
      #print(beta.all[j,87:93])
      BIC[j] <- sum(exp(W%*%beta.all[j,]+offset-as.numeric(t(beta.all[j,])%*%omega%*%beta.all[j,]/2))-Y*(W%*%beta.all[j,])+offset)+sum(beta.all[j,]!=0)*max(log(n),log(log(n))*log(p))
    }
    
    # lambda <- lambdalist[which.min(BIC)]
    # 
    # beta <- SCAD_ADMM_re(X=X, Y=Y, N=N, beta0=beta0, lambda=lambda, rho=rho, a=a, 
    #  
    #                err=err)
    beta <- beta.all[which.min(BIC), ]
  }
  
  return(beta)
}

####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
## Estimate parameters using ADMM for high-dimensional logistic regression 
## without the linear constraints: C beta=b
####################################################################################
## Objective function:
## f(Y,X,beta)+p_a(lambda,theta)+rho v^T (beta-theta)+0.5rho ||beta-theta||_2^2
## f: the loglikelihood
## p: The SCAD penalty function 
## a: parameter in the SCAD penalty, default=3.7
## Y: the response vector
## X: the design matrix
## lambda: the regularization parameter
## rho: regularization parameter in ADMM, we set rho=1
## v: the Langrange multiplier
######################################################

##update_beta
one_step_update_beta <- function(W, Y, offset, omega, beta, theta, v, rho=1){
  n <- length(Y)
  p <- dim(W)[2] 
  
  first_term=crossprod(W/n, Y)
  exp_term1=exp(W%*%beta+offset-as.numeric(t(beta)%*%omega%*%beta/2))
  exp_term2=W-t(matrix(rep(t(beta)%*%omega,n),p,n))
  exp_term=crossprod(exp_term2/n, exp_term1)
  ##p-by-1
  DL1_likelihood=-(first_term-exp_term)
  
  exp_term_DL21=t(exp_term2)%*%diag(c(exp_term1))%*%exp_term2
  exp_term_DL22=sum(exp_term1)*omega
  ##p-by-p
  DL2_likelihood=(exp_term_DL21-exp_term_DL22)/n
  
  beta_new <- (rho*theta-v-DL1_likelihood+  DL2_likelihood%*%beta)
  
  if (p<n){
    #beta_new <- solve( crossprod(X,as.vector(lambda)*X)/n+rho*diag(p), beta_new )
    B0 <- DL2_likelihood+rho*diag(p)
    beta_new1 = ginv(t(B0) %*% B0, tol = mytol) %*% t(B0) %*%  (beta_new)
    #print(beta_new1  )
    #beta_new1= svdB$u %*% diag( svdB$d^(-1)) %*% t(svdB$v) %*% svdB$u %*% diag( svdB$d^(1/2)) %*% t(svdB$v) %*% (beta_new)
    #beta_new1 <- ginv(B0)%*%(beta_new)
    #try(solve(B0, beta_new))
    #if(inherits(beta_new1, "try-error"))
    #  beta_new1 <- ginv(DL2_likelihood+rho*diag(p))%*%beta_new
    beta_new <- beta_new1
  }
  else{
    # using woodbury formula when p is large
    A<-rho*diag(p)-omega*mean(exp_term1)
    U<-t(W-t(matrix(rep(t(beta)%*%omega,n),p,n)))
    C0<-diag(n)/n
    V<-as.vector(exp_term1)*(W-t(matrix(rep(t(beta)%*%omega,n),p,n)))
    invA<-ginv(A,tol=0.0001)
    #try(solve(A))
    #if(inherits(invA, "try-error")){
    #  invA <- ginv(A)
    #}
    invmain<-ginv(n*diag(n)+V%*%invA%*%U,tol=0.0001)
    #try(solve(n*diag(n)+V%*%invA%*%U))
    #if(inherits(invmain, "try-error")){
    #  invmain <- ginv(n*diag(n)+V%*%invA%*%U)
    #}
    beta_new<-(invA-invA%*%U%*%invmain%*%V%*%invA)%*%beta_new
  }
  
  return (beta_new)
}

# ADMM step: update beta
update_beta <- function(W, Y, offset,omega, beta, theta, v, rho=1){
  diff <- 1
  count <- 0
  # estimate beta by Newton Raphson algorithm
#  while(diff>1e-4&&count<=200){
    beta_new <- one_step_update_beta(W, Y,offset, omega, beta, theta, v, rho=rho)
    diff <- sqrt(sum((beta_new-beta)^2))
    
    if (is.na(diff)||is.nan(diff)||is.infinite(diff))
      return(rep(0, length(beta)))
    
    beta <- beta_new
    count <- count+1
 # }
  
  return (beta)
}
#ADMM step: update theta
#SCAD thresholding
# update_theta <- function(beta, v, lambda, rho=1, a=3.7){
#   theta <- v+beta
#   # when absolute value smaller than lambda
#   theta[abs(theta)<=lambda] <- 0
#   # when larger than lambda but smaller than 2lambda
#   index <- abs(theta)>lambda & abs(theta)<=2*lambda
#   theta[index] <- sign(theta[index])*(abs(theta[index])-lambda)
#   # when larger than 2lambda but smaller than a lambda
#   index <- abs(theta)>2*lambda & abs(theta)<=a*lambda
#   theta[index] <- ((a-1)*theta[index]-sign(theta[index])*a*lambda)/(a-2)
# 
#   return (theta)
# }
update_theta <- function(beta, v, N, lambda, rho=1, a=3.7){
  theta <- v+beta
  #N <- c(TRUE, TRUE, rep(FALSE, length(theta)-2))
  # when absolute value smaller than lambda
  #print(theta)
  theta[(!N) & abs(theta)<=lambda] <- 0
  # when larger than lambda but smaller than 2lambda
  index <- (!N) & abs(theta)>lambda & abs(theta)<=2*lambda
  theta[index] <- sign(theta[index])*(abs(theta[index])-lambda)
  # when larger than 2lambda but smaller than a lambda
  index <- (!N) & abs(theta)>2*lambda & abs(theta)<=a*lambda
  theta[index] <- ((a-1)*theta[index]-sign(theta[index])*a*lambda)/(a-2)
  #print(theta)

  return (theta)
}

# calculate the primal and dual residuals
res <- function(beta_new, theta_new, beta_old, theta_old, rho=1){
  # the primal residual
  r <- beta_new-theta_new
  # the dual residual 
  s <- rho*(theta_new-theta_old)
  
  return (list(r=r, s=s))
}

# ADMM algorithm
SCAD_ADMM_uncon <- function(W, Y, offset,N, omega, beta0, lambda, R, rho=1, a=3.7, err=0.5e-3){
		#print(beta0)  
  # initial value of parameters
  beta_old <- beta0
  theta_old <- beta0
  v <- 0
  iter <- 0
  
  eps <- list(r=1, s=1)
  while ((norm(as.matrix(eps$r), "F")>=err || (norm(as.matrix(eps$s), "F")>=err))&&iter<=maxitr)#
  {
    # ADMM update
    # update v
    v <- v+ gamma * (beta_old-theta_old)
    
    # update beta
    beta_new <- update_beta(W=W, Y=Y, offset=offset, omega=omega, beta=beta_old, theta=theta_old, v=v, rho=rho)
    
    # update theta
    
    theta_new <- update_theta(beta=beta_new, v=v, N=N, lambda=lambda, rho=rho, a=a)
    #print(theta_new)
    #project
    theta_new<-euclidean_proj_l1ball(theta_new,s=R*sqrt(2))
    theta_new=shrink(theta_new,R)
    
    # calculate the residuals
    eps <- res(beta_new, theta_new, beta_old, theta_old, rho=rho)
    #print(norm(as.matrix(eps$r), "F"))
    #print(norm(as.matrix(eps$s), "F"))
    #print(beta_new)
    # beta_old and theta_old
    beta_old <- beta_new
    theta_old <- theta_new
    
    iter <- iter+1
    #print(iter)
  }
  
  return(theta_new)
}
##############################################################################
## K-folded cross validation, K=5 by default
##############################################################################
cv.SCAD_ADMM_uncon <- function(W, Y, offset, N, omega, beta0, R, K=5, rho=1, a=3.7, err=0.5e-3, tune="BIC", lambdalist){
  
  # potential tuning parameters
  #lambdalist <- exp(seq(-2.5, 0.5, 0.075))
  #lambdalist <- 5
  lambdan <- length(lambdalist)
  
  # data splitting
  n <- length(Y)
  p <- dim(W)[2]
  beta.all <- matrix(0, lambdan, p)
  
  if (tune=="cv"){
    folds <- split(sample(n, n, replace=FALSE), as.factor(1:K))
    
    # calculate MSE for each folds
    MSE <- rep(0, lambdan)
    for (k in 1:K){
      # testing dataset
      W0 <- W[folds[[k]], ]
      Y0 <- Y[folds[[k]]]
      offset0<-offset[folds[[k]]]
      
      # training dataset
      W1 <- W[-folds[[k]], ]
      Y1 <- Y[-folds[[k]]]
      offset1<-offset[-folds[[k]]]
      
      # est the training dataset
      for (j in 1:lambdan){
        beta0 <- SCAD_ADMM_uncon(W=W1, Y=Y1, offset=offset1,N=N, omega=omega, beta0=beta0, lambda=lambdalist[j], R=R,
                                rho=rho, a=a, err=err)
        MSE[j] <- MSE[j]+sum(exp(W%*%beta0+offset-as.numeric(t(beta0)%*%omega%*%beta0/2))-Y*(W%*%beta0+offset))
      }
    }
    
    # take the minimum lambda
    lambda <- lambdalist[which.min(MSE)]
    beta <- SCAD_ADMM_uncon(W=W, Y=Y, N=N, offset=offset, omega=omega, beta0=beta0, lambda=lambda, R=R, rho=rho, a=a, err=err)
  }
  else{
    BIC <- rep(0, lambdan)
    for (j in 1:lambdan){
      beta0 <- SCAD_ADMM_uncon(W=W, Y=Y, offset=offset, N=N, omega=omega, beta0=beta0, lambda=lambdalist[j], R=R, 
                              rho=rho, a=a, err=err)
      BIC[j] <- sum(exp(W%*%beta0+offset-as.numeric(t(beta0)%*%omega%*%beta0/2))-Y*(W%*%beta0+offset))+sum(beta0!=0)*max(log(n),log(log(n))*log(p))
      beta.all[j,] <- beta0
      #print(beta0[87:93])
    }
    
    # take the minimum lambda
     #lambda <- lambdalist[which.min(BIC)]
     #beta <- SCAD_ADMM_uncon(W=W, Y=Y,N=N, omega=omega, beta0=rep(0,dim(W)[2]), lambda=lambda, R=R, rho=rho, a=a, err=err)
    beta <- beta.all[which.min(BIC), ]
    #print(beta)
  }
  
  return(beta)
}
###construct test statistics
#estimate Q: second derivative of likelihood function
Qest <- function(W, offset,omega, beta){
  n <- dim(W)[1] 
  p <- dim(W)[2] 
  
  exp_term1=exp(W%*%beta+offset-as.numeric(t(beta)%*%omega%*%beta/2))
  exp_term2=W-t(matrix(rep(t(beta)%*%omega,n),p,n))
  exp_term_DL21=t(exp_term2)%*%diag(c(exp_term1))%*%exp_term2
  exp_term_DL22=sum(exp_term1)*omega
  ##p-by-p
  DL2_likelihood=(exp_term_DL21-exp_term_DL22)/n
  return(DL2_likelihood)}

## Qest <- function(W, omega, beta){
##   n <- dim(W)[1] 
##   p <- dim(W)[2] 
  
##   exp_term1=exp(W%*%beta-as.numeric(t(beta)%*%omega%*%beta/2))
##   allW = W- matrix(omega%*%beta, n, p, byrow = T)
##   allW = diag(exp_term1^(1/2)) %*% allW
##   exptermDL21 = (t(allW) %*% allW)/n
##   exp_term_DL22 = mean(exp_term1) * omega
##   ##p-by-p
##   DL2_likelihood=(exp_term_DL21-exp_term_DL22)
##   return(DL2_likelihood)}

#estimate Sigma: variance of score function

Sigest<-function(W, Y, offset, omega, beta){
   first_term=Y*W
   exp_term1=exp(W%*%beta+offset-as.numeric(t(beta)%*%omega%*%beta/2))
   exp_term2=W-t(matrix(rep(t(beta)%*%omega,n),p,n))
   exp_term=as.numeric(exp_term1)*exp_term2
   Sig=t(first_term-exp_term)%*%(first_term-exp_term)/n
   return(Sig)
 }

# Sigest<-function(W, Y, omega, beta){
#   Sig=matrix(0,p,p)
#   for (i in 1:n){
#     Wi=W[i,]
#     Yi=Y[i]
#     term1=Wi%*%beta
#     term2=t(beta)%*%omega%*%beta
#     term=Yi*Wi-as.numeric(exp(term1-term2/2))*(Wi-omega%*%beta)
#     ##term3=(Wi-omega%*%beta)%*%t((Wi-omega%*%beta))
#     Sig_update=term%*%t(term)
#     # print(sigma_update)
#     Sig=Sig+Sig_update
#   }
#   return(Sig/n)
# }
# Sigest1<-function(W, Y, omega, beta){
#   Sig=matrix(0,p,p)
#   for (i in 1:n){
#     Wi=W[i,]
#     Yi=Y[i]
#     term1=Wi%*%beta
#     term2=t(beta)%*%omega%*%beta
#     term3=Wi-omega%*%beta
#     #term=Yi*Wi-as.numeric(exp(term1-term2/2))*(Wi-omega%*%beta)
#     #term3=(Wi-omega%*%beta)%*%t((Wi-omega%*%beta))
#     Sig_update=as.numeric(exp(term1-term2/2))*(term3%*%t(term3)-omega)-
#       as.numeric(exp(2*term1-term2))*(term3%*%t(term3)-omega/2)+
#       as.numeric(exp(term1-term2/2))*omega-as.numeric(exp(2*term1-term2))*omega-
#       as.numeric(exp(2*term1-term2))*omega%*%beta%*%t(term3)-as.numeric(exp(2*term1-term2))*term3%*%t(beta)*omega+
#       as.numeric(exp(2*term1-term2))*term3%*%t(term3)
#     #print(sigma_update)
#     Sig=Sig+Sig_update
#   }
#   return(Sig/n)
# }
Phiest<-function(Sig,Q,beta,C,Index_N){
  index <- union(Index_N, which(beta!=0))
  m <- length(Index_N)
  k <- length(index)-length(Index_N)
  r=dim(C)[1]
  if (r==1){
    Cnew<-t(as.vector(C[,Index_N]))
  }else{Cnew <- as.matrix(C[,Index_N],r,m)}
  Qinv <-  try(solve(Q[index,index]))
  if(inherits(Qinv, "try-error")){
    Qinv <- ginv(Q[index,index], tol = mytol)
  }
  Izero <- cbind(diag(1,m,m),matrix(0,m,k))
  Phi <- Cnew%*%Izero%*%Qinv%*%Sig[index,index]%*%Qinv%*%t(Izero)%*%t(Cnew)
  return(Phi)
}

#first derivative of likelihood function
DL1<-function(W, Y, offset,omega, beta){
first_term=crossprod(W/n, Y)
exp_term1=exp(W%*%beta+offset-as.numeric(t(beta)%*%omega%*%beta/2))
exp_term2=W-t(matrix(rep(t(beta)%*%omega,n),p,n))
exp_term=crossprod(exp_term2/n, exp_term1)
##p-by-1
DL1_likelihood=-(first_term-exp_term)
return(DL1_likelihood)}

