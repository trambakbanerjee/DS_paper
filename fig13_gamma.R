
# effect size estimation
require(ggplot2)
require(gridExtra)
library(ggpubr)
require(Matrix)
library(REBayes)
library(foreach)
library(doParallel)
library(Rcpp)
library(RSpectra)
library(isotone)
library(truncnorm)


##-------------- n varying --------------

exp1_gamma<-function(n,M,reps){
  
  
  eta.grid<-c(0.001,0.0025,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.4,0.5,0.7,1,1.5)
  
  risk.ds_orc<-matrix(0,reps,length(n))
  risk.ds_mcv<-matrix(0,reps,length(n))
  risk.km<-matrix(0,reps,length(n))
  risk.naive<- matrix(0,reps,length(n))
  lam.orc<-matrix(0,reps,length(n))
  lam.mcv<-lam.orc

  for(j in 1:length(n)){
    
    N<- n[j]
    m<-rep(M,N)
    set.seed(j)
    alpha<-sample(c(1,3,5),N,replace = T)
    
    cl <- makeCluster(6)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("truncnorm","isotone","RcppArmadillo","Rcpp","RcppProgress","REBayes","RSpectra","Rfast"))%dopar%{
      
      sourceCpp('dslib.cpp')
      source('ds_lib_final.R')
      set.seed(13*r)
      betaa<- rgamma(N,7,rate = 1)
      X<-lapply(1:N,function(i) gam.rand(m[i],alpha[i],betaa[i],r*i))
      
      x.sum<- sapply(1:N,function(i) sum(X[[i]]))
     
      #Sample splitting
      set.seed(13*r)
      idx<-sample(1:M,0.5*M,replace = F)
      U<- lapply(1:N, function(i) X[[i]][idx])
      u.sum<- sapply(1:N,function(i) sum(U[[i]]))
      uu<- (m*alpha)/u.sum
      V<- lapply(1:N, function(i) X[[i]][-idx])
      v.sum<- sapply(1:N,function(i) sum(V[[i]]))
      vv<- (m*alpha)/v.sum
     
      # Mahalanobis distance
      dat<- cbind(x.sum,alpha)
      dat_u<- cbind(u.sum,alpha)
      dat_v<- cbind(v.sum,alpha)
      Omegainv<- cov(dat)
      Omega<-spdinv(Omegainv)
      Omegainv_u<- cov(dat_u)
      Omega_u<-spdinv(Omegainv_u)
      Omegainv_v<- cov(dat_v)
      Omega_v<-spdinv(Omegainv_v)
      K_md<-matrix(0,N,N)
      K_md_u<-matrix(0,N,N)
      K_md_v<-matrix(0,N,N)
      for(i in 1:N){
        d<-mahala(dat,mu = dat[i,],sigma = Omegainv)
        K_md[,i]<- exp(-0.5*d)
        d<-mahala(dat_u,mu = dat_u[i,],sigma = Omegainv_u)
        K_md_u[,i]<- exp(-0.5*d)
        d<-mahala(dat_v,mu = dat_v[i,],sigma = Omegainv_v)
        K_md_v[,i]<- exp(-0.5*d)
      }
      
      #1. Naive
      naive<- (m*alpha)/x.sum
      risk.naive<- sum((betaa-naive)^2)/N
      
      #2. KM
      out.km<- tryCatch(Gammamix(x.sum, v = 300, shape = m*alpha),
                             error=function(e) NULL)
      
      beta.km<-matrix(0,N,1)
      for(i in 1:N){
        
        num<- sum((1/out.km$x)*dgamma(x.sum[i],m[i]*alpha[i],scale = out.km$x)*out.km$y)
        den<- sum(dgamma(x.sum[i],m[i]*alpha[i],scale = out.km$x)*out.km$y)
        beta.km[i]<-num/den
      }
      
      if(is.null(beta.km)){
        beta.km<-naive
      }
      risk.km_r<- sum((betaa-beta.km)^2)/N
      
      #3. DS Oracle and DS MCV for Mahalanobis Distance
      rho_1<- matrix(0,N,length(eta.grid))
      rho_2<- rho_1
      rho_3<- rho_1
      loss_true_md<-matrix(0,length(eta.grid),1)
      loss_mcv_md<- matrix(0,length(eta.grid),1)
      for(jj in 1:length(eta.grid)){
        
        out_1<-ksd_md_cpp(dat,K_md,Omega[1,1],Omega[1,2],Omega[2,2],eta.grid[jj])
        out_2<-ksd_md_cpp(dat_u,K_md_u,Omega_u[1,1],Omega_u[1,2],Omega_u[2,2],eta.grid[jj])
        out_3<-ksd_md_cpp(dat_v,K_md_v,Omega_v[1,1],Omega_v[1,2],Omega_v[2,2],eta.grid[jj])
        
        rho_1[,jj]<-out_1[,1]
       
        deltaDS_rate<- (m*alpha-1)/x.sum-rho_1[,jj]
        loss_true_md[jj]<-sum((betaa-deltaDS_rate)^2)/N
        
        rho_2[,jj]<-out_2[,1]
        deltaDS_rate_2<- (m*alpha-1)/u.sum-rho_2[,jj]
        rho_3[,jj]<-out_3[,1]
        deltaDS_rate_3<- (m*alpha-1)/v.sum-rho_3[,jj]
        loss_mcv_md[jj]<-0.5*(sum((vv-deltaDS_rate_2)^2)/N+sum((uu-deltaDS_rate_3)^2)/N)
      }
      
      ii<- min(which(loss_true_md==min(loss_true_md,na.rm = T)))
      eta_orc_md<- eta.grid[ii]
      deltaDS_rate<- (m*alpha-1)/x.sum-rho_1[,ii]
      risk.ds_md_orc_r<- sum((betaa-deltaDS_rate)^2)/N
      
      ii<- min(which(loss_mcv_md==min(loss_mcv_md,na.rm = T)))
      eta_mcv_md<- eta.grid[ii]
      out_mcv<-ksd_md_cpp(dat,K_md,Omega[1,1],Omega[1,2],Omega[2,2],
                          ((N<=400)*(1/sqrt(N))+(N>400)*eta_mcv_md))
      rho_mcv<-out_mcv[,1]
      deltaDS_rate<- (m*alpha-1)/x.sum-rho_mcv
      risk.ds_md_mcv_r<- sum((betaa-deltaDS_rate)^2)/N
      
      return(list("risk.naive"=risk.naive,
                  "risk.km"=risk.km_r,
                  "risk.ds_md_orc"=risk.ds_md_orc_r,
                  "risk.ds_md_mcv"=risk.ds_md_mcv_r,
                  "eta_orc_md"=eta_orc_md,
                  "eta_mcv_md"=eta_mcv_md))
    }
    
    stopCluster(cl)
    registerDoSEQ()
    
    risk.ds_orc[,j]<-sapply(1:reps,function(i) result[[i]]$risk.ds_md_orc)
    risk.ds_mcv[,j]<-sapply(1:reps,function(i) result[[i]]$risk.ds_md_mcv)
    risk.naive[,j]<-sapply(1:reps,function(i) result[[i]]$risk.naive)
    risk.km[,j]<-sapply(1:reps,function(i) result[[i]]$risk.km)
    
    lam.orc[,j]<-sapply(1:reps,function(i) result[[i]]$eta_orc_md)
    lam.mcv[,j]<-sapply(1:reps,function(i) result[[i]]$eta_mcv_md)
    print(j)
  }
  result<-list("risk.ds_orc"=risk.ds_orc,"risk.ds_mcv"=risk.ds_mcv,
               "risk.km"=risk.km,"risk.naive"=risk.naive,
               "lam.orc"=lam.orc,
               "lam.mcv"=lam.mcv)
}
exp2_gamma<-function(n,M,reps){
  
  
  eta.grid<-c(0.001,0.0025,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.4,0.5,0.7,1,1.5)
  
  risk.ds_orc<-matrix(0,reps,length(n))
  risk.ds_mcv<-matrix(0,reps,length(n))
  risk.km<-matrix(0,reps,length(n))
  risk.naive<- matrix(0,reps,length(n))
  lam.orc<-matrix(0,reps,length(n))
  lam.mcv<-lam.orc
  
  for(j in 1:length(n)){
    
    N<- n[j]
    m<-rep(M,N)
    set.seed(j)
    alpha<-sample(c(1,2,3),N,replace = T)
    
    cl <- makeCluster(6)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("truncnorm","isotone","RcppArmadillo","Rcpp","RcppProgress","REBayes","RSpectra","Rfast"))%dopar%{
      
      sourceCpp('dslib.cpp')
      source('ds_lib_final.R')
      set.seed(13*r)
      betaa<- rchisq(N,5)
      betaa[betaa<0.1]<- 0.1
      X<-lapply(1:N,function(i) gam.rand(m[i],alpha[i],betaa[i],r*i))
      
      x.sum<- sapply(1:N,function(i) sum(X[[i]]))
      
      #Sample splitting
      set.seed(13*r)
      idx<-sample(1:M,0.5*M,replace = F)
      U<- lapply(1:N, function(i) X[[i]][idx])
      u.sum<- sapply(1:N,function(i) sum(U[[i]]))
      uu<- (m*alpha)/u.sum
      V<- lapply(1:N, function(i) X[[i]][-idx])
      v.sum<- sapply(1:N,function(i) sum(V[[i]]))
      vv<- (m*alpha)/v.sum
      
      # Mahalanobis distance
      dat<- cbind(x.sum,alpha)
      dat_u<- cbind(u.sum,alpha)
      dat_v<- cbind(v.sum,alpha)
      Omegainv<- cov(dat)
      Omega<-spdinv(Omegainv)
      Omegainv_u<- cov(dat_u)
      Omega_u<-spdinv(Omegainv_u)
      Omegainv_v<- cov(dat_v)
      Omega_v<-spdinv(Omegainv_v)
      K_md<-matrix(0,N,N)
      K_md_u<-matrix(0,N,N)
      K_md_v<-matrix(0,N,N)
      for(i in 1:N){
        d<-mahala(dat,mu = dat[i,],sigma = Omegainv)
        K_md[,i]<- exp(-0.5*d)
        d<-mahala(dat_u,mu = dat_u[i,],sigma = Omegainv_u)
        K_md_u[,i]<- exp(-0.5*d)
        d<-mahala(dat_v,mu = dat_v[i,],sigma = Omegainv_v)
        K_md_v[,i]<- exp(-0.5*d)
      }
      
      #1. Naive
      naive<- (m*alpha)/x.sum
      risk.naive<- sum((betaa-naive)^2)/N
      
      #2. KM
      out.km<- tryCatch(Gammamix(x.sum, v = 300, shape = m*alpha),
                        error=function(e) NULL)
      
      beta.km<-matrix(0,N,1)
      for(i in 1:N){
        
        num<- sum((1/out.km$x)*dgamma(x.sum[i],m[i]*alpha[i],scale = out.km$x)*out.km$y)
        den<- sum(dgamma(x.sum[i],m[i]*alpha[i],scale = out.km$x)*out.km$y)
        beta.km[i]<-num/den
      }
      
      if(is.null(beta.km)){
        beta.km<-naive
      }
      risk.km_r<- sum((betaa-beta.km)^2)/N
      
      #3. DS Oracle and DS MCV for Mahalanobis Distance
      rho_1<- matrix(0,N,length(eta.grid))
      rho_2<- rho_1
      rho_3<- rho_1
      loss_true_md<-matrix(0,length(eta.grid),1)
      loss_mcv_md<- matrix(0,length(eta.grid),1)
      for(jj in 1:length(eta.grid)){
        
        out_1<-ksd_md_cpp(dat,K_md,Omega[1,1],Omega[1,2],Omega[2,2],eta.grid[jj])
        out_2<-ksd_md_cpp(dat_u,K_md_u,Omega_u[1,1],Omega_u[1,2],Omega_u[2,2],eta.grid[jj])
        out_3<-ksd_md_cpp(dat_v,K_md_v,Omega_v[1,1],Omega_v[1,2],Omega_v[2,2],eta.grid[jj])
        
        rho_1[,jj]<-out_1[,1]
        
        deltaDS_rate<- (m*alpha-1)/x.sum-rho_1[,jj]
        loss_true_md[jj]<-sum((betaa-deltaDS_rate)^2)/N
        
        rho_2[,jj]<-out_2[,1]
        deltaDS_rate_2<- (m*alpha-1)/u.sum-rho_2[,jj]
        rho_3[,jj]<-out_3[,1]
        deltaDS_rate_3<- (m*alpha-1)/v.sum-rho_3[,jj]
        loss_mcv_md[jj]<-0.5*(sum((vv-deltaDS_rate_2)^2)/N+sum((uu-deltaDS_rate_3)^2)/N)
      }
      
      ii<- min(which(loss_true_md==min(loss_true_md,na.rm = T)))
      eta_orc_md<- eta.grid[ii]
      deltaDS_rate<- (m*alpha-1)/x.sum-rho_1[,ii]
      risk.ds_md_orc_r<- sum((betaa-deltaDS_rate)^2)/N
      
      ii<- min(which(loss_mcv_md==min(loss_mcv_md,na.rm = T)))
      eta_mcv_md<- eta.grid[ii]
      out_mcv<-ksd_md_cpp(dat,K_md,Omega[1,1],Omega[1,2],Omega[2,2],
                          eta_mcv_md)
      rho_mcv<-out_mcv[,1]
      deltaDS_rate<- (m*alpha-1)/x.sum-rho_mcv
      risk.ds_md_mcv_r<- sum((betaa-deltaDS_rate)^2)/N
      
      return(list("risk.naive"=risk.naive,
                  "risk.km"=risk.km_r,
                  "risk.ds_md_orc"=risk.ds_md_orc_r,
                  "risk.ds_md_mcv"=risk.ds_md_mcv_r,
                  "eta_orc_md"=eta_orc_md,
                  "eta_mcv_md"=eta_mcv_md))
    }
    
    stopCluster(cl)
    registerDoSEQ()
    
    risk.ds_orc[,j]<-sapply(1:reps,function(i) result[[i]]$risk.ds_md_orc)
    risk.ds_mcv[,j]<-sapply(1:reps,function(i) result[[i]]$risk.ds_md_mcv)
    risk.naive[,j]<-sapply(1:reps,function(i) result[[i]]$risk.naive)
    risk.km[,j]<-sapply(1:reps,function(i) result[[i]]$risk.km)
    
    lam.orc[,j]<-sapply(1:reps,function(i) result[[i]]$eta_orc_md)
    lam.mcv[,j]<-sapply(1:reps,function(i) result[[i]]$eta_mcv_md)
    print(j)
  }
  result<-list("risk.ds_orc"=risk.ds_orc,"risk.ds_mcv"=risk.ds_mcv,
               "risk.km"=risk.km,"risk.naive"=risk.naive,
               "lam.orc"=lam.orc,
               "lam.mcv"=lam.mcv)
}


#--------------------------------------------------------------------

n<-seq(100,1000,length.out = 10)
reps<-200
M<- 30

out1_m30<-exp1_gamma(n,M,reps)
out2_m30<-exp2_gamma(n,M,reps)


plotdata1<- as.data.frame(c(colMeans(out1_m30$risk.km)/colMeans(out1_m30$risk.ds_orc),
                            colMeans(out1_m30$risk.ds_mcv)/colMeans(out1_m30$risk.ds_orc),
                            colMeans(out1_m30$risk.naive)/colMeans(out1_m30$risk.ds_orc)))
names(plotdata1)<-"Risk"
plotdata1$type<-as.factor(rep(c('NPMLE',
                                'DS','MLE'),each=length(n)))
plotdata1$n<- rep(n,3)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=n,y=Risk,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=n,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+scale_x_continuous(breaks = n)+ylab("Relative Risk")+
  #scale_y_continuous(breaks = c(0.85,0.9,0.95,1.0,1.05,1.1))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata2<- as.data.frame(c(colMeans(out2_m30$risk.km)/colMeans(out2_m30$risk.ds_orc),
                            colMeans(out2_m30$risk.ds_mcv)/colMeans(out2_m30$risk.ds_orc),
                            colMeans(out2_m30$risk.naive)/colMeans(out2_m30$risk.ds_orc)))
names(plotdata2)<-"Risk"
plotdata2$type<-as.factor(rep(c('NPMLE',
                                'DS','MLE'),each=length(n)))
plotdata2$n<- rep(n,3)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=n,y=Risk,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=n,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+scale_x_continuous(breaks = n)+ylab("Relative Risk")+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))



ggarrange(g1,g2,ncol=2,nrow=1,common.legend = TRUE)


save.image(paste(getwd(),'fig13_gamma.RData',sep=''))
