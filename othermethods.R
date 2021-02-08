

library(MASS)
library(isotone)

#-------------------------------------------------------------
## group-linear and XKB functions

## spherically symmetric estimator with c_n = c^*_n
spher <- function(x.,v.){
  n. <- length(x.)
  if ( (n.==1) | (var(x.)==0) ) x. else {
    cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)
    bhat <- min( cstar*mean(v.)/var(x.), 1 )
    x. - bhat*(x. - mean(x.))
  }
}


## spherically symmetric estimator with c_n = c^*_n, shrinkage toward zero
spher.zero <- function(x.,v.){
  n. <- length(x.)
  cstar <- max( 1-2*( max(v.)/mean(v.) )/n., 0)
  bhat <- min( cstar*mean(v.)/mean(x.^2), 1 )
  (1- bhat)*x.
}

## function that returns the common bhat (replicated)
spher.bhat <- function(x.,v.){
  n. <- length(x.)
  if ( (n.==1) | (var(x.)==0) ) x. else {
    cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)
    bhat <- min( cstar*mean(v.)/var(x.), 1 )
    return(rep(bhat,n.))
  }
}

## group-linear estimator

grouplinear <- function( x,v,nbreak=floor(length(x)^(1/3)) ){  # default: bin log(v) into same NUMBER (=n^(1/3) of intervals
  n <- length(x)
  splitby=cut(log(v),breaks=nbreak, labels=F)
  xsub <- split(x,splitby)
  vsub <- split(v,splitby)
  indexsub <- split(1:n,splitby)
  thetahatsub <- mapply(spher,xsub,vsub)
  indexsub.unlist <- as.vector( unlist(indexsub) )
  thetahatsub.unlist <- as.vector( unlist(thetahatsub) )
  thetahat <- thetahatsub.unlist[order(indexsub.unlist)]	
  return(thetahat)
}

grouplinear.zero <- function( x,v,nbreak=floor(length(x)^(1/3)) ){  # default: bin log(v) into same NUMBER (=n^(1/3) of intervals
  n <- length(x)
  splitby=cut(log(v),breaks=nbreak, labels=F)
  xsub <- split(x,splitby)
  vsub <- split(v,splitby)
  indexsub <- split(1:n,splitby)
  thetahatsub <- mapply(spher.zero,xsub,vsub)
  indexsub.unlist <- as.vector( unlist(indexsub) )
  thetahatsub.unlist <- as.vector( unlist(thetahatsub) )
  thetahat <- thetahatsub.unlist[order(indexsub.unlist)]	
  return(thetahat)
}

thetahat.M <- function(X,A){
  lambda.sure <- ifelse( g(0,X=X,A=A)*g(max(A)*1000,X=X,A=A) < 0, uniroot(g,c(0,max(A)*1000),X=X,A=A, tol=1e-9)$root, optim(c(mean(   pmax(  ( X-mean(X) )^2 - 	A,0  )   ), mean(X)),f,X=X,A=A, method = "L-BFGS-B",lower=c(0,-Inf))$par[1] )
  mu.sure <- sum( A^2/(lambda.sure+A)^2 * X ) / sum( A^2/(lambda.sure+A)^2 )
  thetahat(X,A,lambda.sure,mu.sure)
}

thetahat.G <- function(X,A){
  lambda <- optimize(f.G,lower=0,upper=1000,X=X,A=A)$minimum
  thetahat(X,A,lambda,mean(X))
}


thetahat.SG <- function(X,A){
  p <- length(X)
  fit <- gpava( z = A, y = A * (1-1/p) / (X-mean(X))^2, weights = (X-mean(X))^2, solver = weighted.mean, ties="primary" )
  bhat <- pmin(  pmax( fit$x,0 ),1  )
  (1-bhat) * X + bhat * mean(X)
}

f.G <- function(lambda,X,A){ 
  p <- length(X)
  sum(  ( A/(A+lambda) )^2 * (X - mean(X))^2 + A/(A+lambda) * (lambda - A + 2/p * A)  )
}

thetahat <- function(X,A,lambda,mu){
  lambda/(lambda+A) * X + A/(lambda+A) * mu
}

# d/dlambda{SURE(lambda,mu=mu.hat.SURE(lambda))} (proportional to)
g <- function(lambda,X,A){  
  sum( A^2/(lambda+A)^3 * (X-(  sum( A^2/(lambda+A)^2 * X ) / sum( A^2/(lambda+A)^2 )  ))^2 - A^2/(lambda+A)^2 )
  #equivalent to the following(which is just easier to read):
  #mu <- sum( A^2/(lambda+A)^2 * X ) / sum( A^2/(lambda+A)^2 )
  #sum( A^2/(lambda+A)^3 * (X-mu)^2 - A^2/(lambda+A)^2 )
}

#------------------------------------------------

# EBMN-moderated -t (Stephens 2019)

ebnm_mod_t<-function(v){
  
  # Use lmFit to estimate effect sizes betahat
  lim = lmFit(v)
  # Then use eBayes to shrink the standard errors
  lim = eBayes(lim)
  # Get betahat, the shrunk s.e. and moderated d.f
  betahat = c(lim$coefficients)
  se = c(sqrt(lim$s2.post)) # EB shrunk s.e.#c(lim$stdev.unscaled*sqrt(lim$s2.post)) # EB shrunk s.e.
  df = lim$df.total[1]
  # Fit ash with the shrunk s.e. and moderated d.f
  fit = ash(betahat, se, df=df)
  return(list('thetahat'=fit$result$betahat,'sigmahat'=se))
}
