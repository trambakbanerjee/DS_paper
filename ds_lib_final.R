
library(MASS)
library(Rfast)

mvnorm.rand<-function(nobs,mean,sigma,setseed){
  
  set.seed(setseed)
  return(mvrnorm(nobs,mu=mean,Sigma=sigma))
}

norm.rand<-function(nobs,mean,sd,setseed){
  
  set.seed(setseed)
  return(rnorm(nobs,mean,sd))
}
lognorm.rand<-function(nobs,mean,sd,setseed){
  
  set.seed(setseed)
  return(exp(rnorm(nobs,mean,sd)))
}
gam.rand<-function(nobs,aa,bb,setseed){
  
  set.seed(setseed)
  return(rgamma(nobs,aa,rate=bb))
}
weibull.rand<-function(nobs,aa,bb,setseed){
  
  set.seed(setseed)
  return(rweibull(nobs,aa,scale=bb))
}

unif.rand<-function(nobs,aa,bb,setseed){
  
  set.seed(setseed)
  return(runif(nobs,aa,bb))
}
normexp.rand<-function(nobs,aa,mean,sd,setseed){
  
  set.seed(setseed)
  u<- rexp(nobs,rate=aa)
  set.seed(10*setseed)
  v<-rnorm(nobs,mean,sd)
  return(u+v)
}
normt.rand<-function(nobs,dd,mean,sd,setseed){
  
  set.seed(setseed)
  u<- rt(nobs,dd,0)
  set.seed(10*setseed)
  v<-rnorm(nobs,mean,sd)
  return(u+v)
}
normlap.rand<-function(nobs,dd,mean,sd,setseed){
  
  set.seed(setseed)
  u<- rlaplace(nobs,0,dd)
  set.seed(10*setseed)
  v<-rnorm(nobs,mean,sd)
  return(u+v)
}
t.rand<-function(nobs,dd,setseed){
  
  set.seed(setseed)
  u<- rt(nobs,dd,0)
  return(u)
}
beta.rand<-function(nobs,aa,bb,setseed){
  
  set.seed(setseed)
  return(rbeta(nobs,aa,bb))
}
lap.rand<-function(nobs,dd,setseed){
  
  set.seed(setseed)
  u<- rlaplace(nobs,0,dd)
  return(u)
}
cauchy.rand<-function(nobs,setseed){
  
  set.seed(setseed)
  u<- rcauchy(nobs,0,0.005)
  return(u)
}
pois.rand<-function(nobs,ll,setseed){
  
  set.seed(setseed)
  u<- rpois(nobs,ll)
  return(u)
}

