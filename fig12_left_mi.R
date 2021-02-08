
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



##-------------- Fixed n and sigma varying --------------

exp_5<-function(n,sigma2_u,reps){
  
  
  alpha<-0.5
  eta.grid<-c(4:10)
  
  risk.ds_orc<-matrix(0,reps,length(sigma2_u))
  risk.ds_mcv<-matrix(0,reps,length(sigma2_u))
  risk.gl<-matrix(0,reps,length(sigma2_u))
  risk.xkb.g<-matrix(0,reps,length(sigma2_u))
  risk.xkb.sg<-matrix(0,reps,length(sigma2_u))
  risk.xkb.m<-matrix(0,reps,length(sigma2_u))
  risk.km<-matrix(0,reps,length(sigma2_u))
  lam.orc<-matrix(0,reps,length(sigma2_u))
  lam.mcv<-lam.orc
 
  N<- n
  set.seed(1)
  m<-sample(c(10,15,20,25,30,35,40,45,50),N,replace = T)
  
  for(j in 1:length(sigma2_u)){
    
    uu<-sigma2_u[j]
    
    cl <- makeCluster(6)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("truncnorm","isotone","RcppArmadillo","Rcpp","RcppProgress","REBayes","RSpectra","Rfast"))%dopar%{
      
      sourceCpp('dslib.cpp')
      source('ds_lib_final.R')
      source('othermethods.R')
      set.seed(10*r^2)
      p<-runif(N)
      set.seed(r)
      A<-rtruncnorm(N,a=0.01,mean = uu,sd = 1)
      set.seed(100*r)
      theta<- (p<=0.8)*c(rmvnorm(1,A/4,0.25*diag(N)))+(p>0.8)*c(rmvnorm(1,A,diag(N)))
      X<-lapply(1:N,function(i) unif.rand(m[i],theta[i]-sqrt(3*A[i]),theta[i]+sqrt(3*A[i]),r*i))
      
      xbar<- sapply(1:N,function(i) mean(X[[i]]))
      S2.sample<- sapply(1:N,function(i)var(X[[i]]))
      S.sample<-sqrt(S2.sample)
      
      #Sample splitting
      eps<-lapply(1:N,function(i) norm.rand(m[i],0,S.sample[i],10*r*i))
      ubar<-sapply(1:N,function(i) mean(X[[i]]+alpha*eps[[i]]))
      ebar<-sapply(1:N,function(i) mean(eps[[i]]))
      S2.u<- sapply(1:N,function(i)var(X[[i]]+alpha*eps[[i]]))
      vbar<-sapply(1:N,function(i) mean(X[[i]]-(1/alpha)*eps[[i]]))
      S2.v<- sapply(1:N,function(i)var(X[[i]]-(1/alpha)*eps[[i]]))
      
      # Mahalanobis distance
      dat<- cbind(xbar,S2.sample,m)
      dat_u<- cbind(ubar,S2.u,m)
      Omegainv<- cov(dat)
      Omega<-spdinv(Omegainv)
      Omegainv_u<- cov(dat_u)
      Omega_u<-spdinv(Omegainv_u)
      K_md<-matrix(0,N,N)
      K_md_u<-matrix(0,N,N)
      for(i in 1:N){
        d<-mahala(dat,mu = dat[i,],sigma = Omegainv)
        K_md[,i]<- exp(-0.5*d)
        d<-mahala(dat_u,mu = dat_u[i,],sigma = Omegainv_u)
        K_md_u[,i]<- exp(-0.5*d)
      }
      
      #1. Group Linear
      thetahat.gl<- grouplinear(xbar,S2.sample/m)
      risk.gl_r<- sum((theta-thetahat.gl)^2)/N
      
      #2. XKB - 3 versions
      thetahat.xkb.g<- thetahat.G(xbar,S2.sample/m)
      risk.xkb_g<- sum((theta-thetahat.xkb.g)^2)/N
      
      thetahat.xkb.sg<- thetahat.SG(xbar,S2.sample/m)
      risk.xkb_sg<- sum((theta-thetahat.xkb.sg)^2)/N
      
      thetahat.xkb.m<- thetahat.M(xbar,S2.sample/m)
      risk.xkb_m<- sum((theta-thetahat.xkb.m)^2)/N
      
      #3. KM
      X.vec<-unlist(X)
      id.km<-rep(1:N,m)
      thetahat.km<- tryCatch(WGLVmix(X.vec,id.km,rep(1,length(X.vec)))$dx,
                             error=function(e) NULL)
      if(is.null(thetahat.km)){
        thetahat.km<-xbar
      }
      risk.km_r<- sum((theta-thetahat.km)^2)/N
      
      #4. DS Oracle and DS MCV for Mahalanobis Distance
      rho1_2<- matrix(0,N,length(eta.grid))
      rho2_2<- rho1_2
      rho1_4<- rho1_2
      rho2_4<- rho2_2
      loss_mcv_md<-matrix(0,length(eta.grid),1)
      loss_true_md<-matrix(0,length(eta.grid),1)
      for(jj in 1:length(eta.grid)){
        
        out_2<-ksd_3d_md_cpp(dat,K_md,Omega[1,1],Omega[1,2],Omega[1,3],
                             Omega[2,2],Omega[2,3],eta.grid[jj])
        out_4<-ksd_3d_md_cpp(dat_u,K_md_u,Omega_u[1,1],Omega_u[1,2],Omega_u[1,3],
                             Omega_u[2,2],Omega_u[2,3],eta.grid[jj])
        
        rho1_2[,jj]<-out_2[,1]
        rho2_2[,jj]<-out_2[,2]
        
        tauhat<- (m-3)/((m-1)*S2.sample)-(2/(m-1))*rho2_2[,jj]
        deltaDS<- xbar+(1/(m*tauhat))*rho1_2[,jj]
        loss_true_md[jj]<-sum((theta-deltaDS)^2)/N
        
        rho1_4[,jj]<-out_4[,1]
        rho2_4[,jj]<-out_4[,2]
        
        tauhat<- (m-3)/((m-1)*S2.u)-(2/(m-1))*rho2_4[,jj]
        deltaDS<- ubar+(1/(m*tauhat))*rho1_4[,jj]
        loss_mcv_md[jj]<-sum((vbar-deltaDS)^2)/N
      }
      
      ii<- min(which(loss_true_md==min(loss_true_md,na.rm = T)))
      eta_orc_md<- eta.grid[ii]
      tauhat<- (m-3)/((m-1)*S2.sample)-(2/(m-1))*rho2_2[,ii]
      deltaDS<- xbar+(1/(m*tauhat))*rho1_2[,ii]
      risk.ds_md_orc_r<- sum((theta-deltaDS)^2)/N
      
      ii<- min(which(loss_mcv_md==min(loss_mcv_md,na.rm = T)))
      eta_mcv_md<- eta.grid[ii]
      tauhat<- (m-3)/((m-1)*S2.sample)-(2/(m-1))*rho2_2[,ii]
      deltaDS<- xbar+(1/(m*tauhat))*rho1_2[,ii]
      risk.ds_md_mcv_r<- sum((theta-deltaDS)^2)/N
      
   
      return(list("risk.gl"=risk.gl_r,
                  "risk.xkb.g"=risk.xkb_g,"risk.xkb.sg"=risk.xkb_sg,"risk.xkb.m"=risk.xkb_m,
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
    risk.gl[,j]<-sapply(1:reps,function(i) result[[i]]$risk.gl)
    risk.xkb.g[,j]<-sapply(1:reps,function(i) result[[i]]$risk.xkb.g)
    risk.xkb.sg[,j]<-sapply(1:reps,function(i) result[[i]]$risk.xkb.sg)
    risk.xkb.m[,j]<-sapply(1:reps,function(i) result[[i]]$risk.xkb.m)
    risk.km[,j]<-sapply(1:reps,function(i) result[[i]]$risk.km)
    
    lam.orc[,j]<-sapply(1:reps,function(i) result[[i]]$eta_orc_md)
    lam.mcv[,j]<-sapply(1:reps,function(i) result[[i]]$eta_mcv_md)
    print(j)
  }
  result<-list("risk.ds_orc"=risk.ds_orc,"risk.ds_mcv"=risk.ds_mcv,
               "risk.gl"=risk.gl,"risk.km"=risk.km,"risk.xkbg"=risk.xkb.g,
               "risk.xkbsg"=risk.xkb.sg,"risk.xkbm"=risk.xkb.m,"lam.orc"=lam.orc,
               "lam.mcv"=lam.mcv)
}

#--------------------------------------------------------------------

n<-1000
reps<-50
sigma2_u<- c(0.5,1,1.5,2,2.5,3)

out_mi<-exp_5(n,sigma2_u,reps)


plotdata1<- as.data.frame(c(colMeans(out_mi$risk.gl)/colMeans(out_mi$risk.ds_orc),
                            colMeans(out_mi$risk.km)/colMeans(out_mi$risk.ds_orc),
                            colMeans(out_mi$risk.xkbsg)/colMeans(out_mi$risk.ds_orc),
                            colMeans(out_mi$risk.xkbm)/colMeans(out_mi$risk.ds_orc),
                            colMeans(out_mi$risk.ds_mcv)/colMeans(out_mi$risk.ds_orc)))
names(plotdata1)<-"Risk"
plotdata1$type<-as.factor(rep(c('Grp Linear','NPMLE','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata1$u<- rep(sigma2_u,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

save.image(paste(getwd(),'fig12_left_mi.RData',sep=''))
