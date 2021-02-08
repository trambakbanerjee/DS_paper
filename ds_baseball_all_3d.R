
require(ggplot2)
require(gridExtra)
library(ggpubr)
require(Matrix)
library(REBayes)
library(foreach)
library(doParallel)
library(Rcpp)
library(RSpectra)
library(Rfast)
library(isotone)
library(reshape2)
library(zoo)


sourceCpp('dslib.cpp')
source('ds_lib_final.R')
source('othermethods.R')

# Process data ---------------------------

baseball<-bball
dat.1<-baseball[baseball$year<2012,]
cc<-as.data.frame(table(dat.1$id))
dd<- as.data.frame(cc$Var1[cc$Freq>=5])
names(dd)<-'id'

idx<-which(baseball$id %in% dd$id)
dat.bball<-baseball[idx,]
dat.1<-dat.bball[dat.bball$year<2012,]
players.1<-as.data.frame(unique(as.character(dat.1$name)))
names(players.1)<-'name'

n<-length(players.1$name)
X<-list()
v<-list()
m<-matrix(0,n,1)
S2.sample<-matrix(0,n,1)
for(i in 1:n){
  ii<- as.character(players.1$name[i])
  temp<-dat.1[as.character(dat.1$name) == ii,]
  X[[i]]<- temp$HA
  m[i]<-length(X[[i]])
  v[[i]]<-4*temp$AB
}
xbar<- sapply(1:n,function(i) sum(X[[i]]*v[[i]])/(sum(v[[i]])))
multi<- sapply(1:n,function(i) sum(v[[i]]))
S2.sample<- sapply(1:n,function(i) 
  (1/(m[i]-1))*(sum(v[[i]]*X[[i]]^2)-(xbar[i]^2)*sum(v[[i]])))
S.sample<-sqrt(S2.sample)

dat.2<-dat.bball[dat.bball$year==2012,]
mm <- melt(dat.2, measure.vars = c("AB","H"))
dout <- dcast(mm, name + id+year ~ variable, fun.aggregate = sum, margins = "year")
dout<-dout[dout$year!="(all)",]
theta<- sapply(1:length(unique(dout$name)),
               function(i) asin(sqrt((dout$H[i]+0.25)/(dout$AB[i]+0.5))))
pi<-sapply(1:length(unique(dout$name)),
           function(i) dout$H[i]/dout$AB[i])
idy<-which(players.1$name %in% dout$name)


#------------- Methods----------------------
#1. REBayes with independent priors
X.vec<-unlist(X)
id.km<-dat.1$id
w<-unlist(v)
thetahat.km_ind<- tryCatch(WTLVmix(X.vec,id.km,w)$dx,
                       error=function(e) NULL)
if(is.null(thetahat.km_ind)){
  thetahat.km_ind<-xbar
}
p.km_ind<-(sin(thetahat.km_ind))^2
TSE.km_ind<- sum(((theta-thetahat.km_ind[idy])^2-(0.25/dout$AB))^1)
TSEp.km_ind<- sum((pi-p.km_ind[idy])^2-pi*(1-pi)/dout$AB)
NSE.km_ind<- sum(4*dout$AB*(theta-thetahat.km_ind[idy])^2)

#2. REBayes with dependent priors
thetahat.km_dep<- tryCatch(WGLVmix(X.vec,id.km,w,u = 100, v = 100)$dx,
                           error=function(e) NULL)
if(is.null(thetahat.km_dep)){
  thetahat.km_dep<-xbar
}
p.km_dep<-(sin(thetahat.km_dep))^2
TSE.km_dep<- sum(((theta-thetahat.km_dep[idy])^2-(0.25/dout$AB))^1)
TSEp.km_dep<- sum((pi-p.km_dep[idy])^2-pi*(1-pi)/dout$AB)
NSE.km_dep<- sum(4*dout$AB*(theta-thetahat.km_dep[idy])^2)

#3. Naive using only 2011
loss1.naive<-matrix(0,length(idy),1)
loss2.naive<-matrix(0,length(idy),1)
loss3.naive<-matrix(0,length(idy),1)

for(i in 1:length(idy)){
  
  dat.naive<-dat.1[dat.1$name==dout$name[i],]
  dat.naive<-dat.naive[dat.naive$year==max(dat.naive$year),]
  if(dim(dat.naive)[1]>0){
  mm <- melt(dat.naive, measure.vars = c("AB","H"))
  temp <- dcast(mm, name + id+year ~ variable, fun.aggregate = sum, margins = "year")
  temp<-temp[temp$year!="(all)",]
 
  loss1.naive[i]<- (theta[i]-asin(sqrt((temp$H+0.25)/(temp$AB+0.5))))^2-(0.25/dout$AB[i])
  loss2.naive[i]<- 4*dout$AB[i]*(theta[i]-asin(sqrt((temp$H+0.25)/(temp$AB+0.5))))^2
  loss3.naive[i]<-(pi[i]-(temp$H+0.25)/(temp$AB+0.5))^2-pi[i]*(1-pi[i])/dout$AB[i]
  }
}
TSE.naive<- sum(loss1.naive)
TSEp.naive<- sum(loss3.naive)
NSE.naive<- sum(loss2.naive)

#4. Naive using average of 2011
dat.lag<-dat.1[dat.1$year==2011,]
lag2011<-mean(dat.lag$HA)
plag2011<-(sin(lag2011))^2
TSE.lag<- sum(((theta-lag2011)^2-(0.25/dout$AB))^1)
TSEp.lag<- sum((pi-plag2011)^2-pi*(1-pi)/dout$AB)
NSE.lag<- sum(4*dout$AB*(theta-lag2011)^2)

#5. Naive using weighted historical average (2002-2011)
phist<-(sin(xbar))^2
TSE.hist<- sum(((theta-xbar[idy])^2-(0.25/dout$AB))^1)
TSEp.hist<- sum((pi-phist[idy])^2-pi*(1-pi)/dout$AB)
NSE.hist<- sum(4*dout$AB*(theta-xbar[idy])^2)

#6. Naive using unweighted historical average (2002-2011)
xbar.unwt<- sapply(1:n,function(i) mean(X[[i]]))
phist.unwt<-(sin(xbar.unwt))^2
TSE.histunwt<- sum(((theta-xbar.unwt[idy])^2-(0.25/dout$AB))^1)
TSEp.histunwt<- sum((pi-phist.unwt[idy])^2-pi*(1-pi)/dout$AB)
NSE.histunwt<- sum(4*dout$AB*(theta-xbar.unwt[idy])^2)

#--------------------------------
#1. Group Linear
thetahat.gl<- grouplinear(xbar,S2.sample/multi)
p.gl<-(sin(thetahat.gl))^2
TSE.gl<- sum(((theta-thetahat.gl[idy])^2-(0.25/dout$AB))^2)
TSEp.gl<- sum((pi-p.gl[idy])^2-pi*(1-pi)/dout$AB)
NSE.gl<- sum(4*dout$AB*(theta-thetahat.gl[idy])^2)

#2. XKB - 3 versions
thetahat.xkb.g<- thetahat.G(xbar,S2.sample/multi)
p.xkbg<-(sin(thetahat.xkb.g))^2
TSE.xkbg<- sum(((theta-thetahat.xkb.g[idy])^2-(0.25/dout$AB))^2)
TSEp.xkbg<- sum((pi-p.xkbg[idy])^2-pi*(1-pi)/dout$AB)
NSE.xkbg<- sum(4*dout$AB*(theta-thetahat.xkb.g[idy])^2)

thetahat.xkb.sg<- thetahat.SG(xbar,S2.sample/multi)
p.xkbsg<-(sin(thetahat.xkb.sg))^2
TSE.xkbsg<- sum(((theta-thetahat.xkb.sg[idy])^2-(0.25/dout$AB))^2)
TSEp.xkbsg<- sum((pi-p.xkbsg[idy])^2-pi*(1-pi)/dout$AB)
NSE.xkbsg<- sum(4*dout$AB*(theta-thetahat.xkb.sg[idy])^2)

thetahat.xkb.m<- thetahat.M(xbar,S2.sample/multi)
p.xkbm<-(sin(thetahat.xkb.m))^2
TSE.xkbm<- sum(((theta-thetahat.xkb.m[idy])^2-(0.25/dout$AB))^2)
TSEp.xkbm<- sum((pi-p.xkbm[idy])^2-pi*(1-pi)/dout$AB)
NSE.xkbm<- sum(4*dout$AB*(theta-thetahat.xkb.m[idy])^2)
#----------------------------------------------

#7. Double Shrinkage  
eta.grid<-seq(0.1,0.5,length.out = 50)

#Sample splitting
alpha<-0.5

eps<-lapply(1:n,function(i) mvnorm.rand(1,rep(0,m[i]),diag(S2.sample[i]/v[[i]]),10*i))
ebar<- sapply(1:n,function(i) sum(eps[[i]]*v[[i]])/(sum(v[[i]])))
ubar<-xbar+alpha*ebar
S2.u<- S2.sample+(alpha^2)*sapply(1:n,function(i) 
  (1/(m[i]-1))*(sum(v[[i]]*eps[[i]]^2)-(ebar[i]^2)*sum(v[[i]])))
vbar<-xbar-(1/alpha)*ebar
S2.v<- S2.sample+(1/alpha^2)*sapply(1:n,function(i) 
  (1/(m[i]-1))*(sum(v[[i]]*eps[[i]]^2)-(ebar[i]^2)*sum(v[[i]])))

# Mahalanobis distance
dat<- cbind(scale(xbar),scale(S2.sample),scale(m))
dat_u<- cbind(scale(ubar),scale(S2.u),scale(m))
Omegainv<- cov(dat)
Omega<-spdinv(Omegainv)
Omegainv_u<- cov(dat_u)
Omega_u<-spdinv(Omegainv_u)
K_md<-matrix(0,n,n)
K_md_u<-matrix(0,n,n)
for(i in 1:n){
  d<-mahala(dat,mu = dat[i,],sigma = Omegainv)
  K_md[,i]<- exp(-0.5*d)
  d<-mahala(dat_u,mu = dat_u[i,],sigma = Omegainv_u)
  K_md_u[,i]<- exp(-0.5*d)
}
#4. DS Oracle (using squared error loss) and DS MCV for Mahalanobis Distance

rho1_2<- matrix(0,n,length(eta.grid))
rho2_2<- rho1_2
rho1_4<- rho1_2
rho2_4<- rho2_2
loss_true_md<-matrix(0,length(eta.grid),1)
loss_mcv_md<- matrix(0,length(eta.grid),1)

for(jj in 1:length(eta.grid)){

    out_2<-ksd_3d_md_cpp(dat,K_md,Omega[1,1],Omega[1,2],Omega[1,3],
                         Omega[2,2],Omega[2,3],eta.grid[jj])
    out_4<-ksd_3d_md_cpp(dat_u,K_md_u,Omega_u[1,1],Omega_u[1,2],Omega_u[1,3],
                         Omega_u[2,2],Omega_u[2,3],eta.grid[jj])

    rho1_2[,jj]<-out_2[,1]
    rho2_2[,jj]<-out_2[,2]

    tauhat<- (m-3)/((m-1)*S2.sample)-(2/(m-1))*rho2_2[,jj]
    deltaDS<- xbar+(1/(multi*tauhat))*rho1_2[,jj]
    loss_true_md[jj]<-sum((theta-deltaDS[idy])^2)

    rho1_4[,jj]<-out_4[,1]
    rho2_4[,jj]<-out_4[,2]

    tauhat<- (m-3)/((m-1)*S2.u)-(2/(m-1))*rho2_4[,jj]
    deltaDS<- ubar+(1/(multi*tauhat))*rho1_4[,jj]
    loss_mcv_md[jj]<-sum((vbar-deltaDS)^2)

  }
  
  ii<- min(which(loss_true_md==min(loss_true_md,na.rm = T)))
  eta_orc_md<- eta.grid[ii]
  
  tauhat<-((m-3))/((m-1)*S2.sample)-(2/(m-1))*rho2_2[,ii]
  tauhat[tauhat<0]<-(m[tauhat<0]-3)/((m[tauhat<0]-1)*S2.sample[tauhat<0])
  tauhat_orc<- tauhat
  deltaDS_orc<- xbar+(1/(multi*tauhat_orc))*rho1_2[,ii]
  phat_orc<-(sin(deltaDS_orc))^2
  
  TSE.ds_md_orc<- sum((theta-deltaDS_orc[idy])^2-(0.25/dout$AB))
  TSEp.ds_md_orc<- sum((pi-phat_orc[idy])^2-pi*(1-pi)/dout$AB)
  NSE.ds_md_orc<- sum(4*dout$AB*(theta-deltaDS_orc[idy])^2)
  
  
  ii<- min(which(loss_mcv_md==min(loss_mcv_md,na.rm = T)))
  eta_mcv_md<- eta.grid[ii]
  tauhat<- (m-3)/((m-1)*S2.sample)-(2/(m-1))*rho2_2[,ii]
  tauhat[tauhat<0]<-(m[tauhat<0]-3)/((m[tauhat<0]-1)*S2.sample[tauhat<0])
  tauhatDS<- tauhat
  deltaDS<- xbar+(1/(multi*tauhatDS))*rho1_2[,ii]
  phat_DS<-(sin(deltaDS))^2
  
  TSE.ds_md_mcv<- sum((theta-deltaDS[idy])^2-(0.25/dout$AB))
  TSEp.ds_md_mcv<- sum((pi-phat_DS[idy])^2-pi*(1-pi)/dout$AB)
  NSE.ds_md_mcv<- sum(4*dout$AB*(theta-deltaDS[idy])^2)
  
  final_table<- as.data.frame(cbind(c(TSE.km_ind,TSE.km_dep,TSE.hist,TSE.histunwt,TSE.lag,TSE.naive,
                                      TSE.ds_md_mcv,TSE.ds_md_orc)/TSE.naive,
                              c(NSE.km_ind,NSE.km_dep,NSE.hist,NSE.histunwt,NSE.lag,NSE.naive,
                                NSE.ds_md_mcv,NSE.ds_md_orc)/NSE.naive,
  c(TSEp.km_ind,TSEp.km_dep,TSEp.hist,TSEp.histunwt,TSEp.lag,TSEp.naive,
    TSEp.ds_md_mcv,TSEp.ds_md_orc)/TSEp.naive))
  
  colnames(final_table)<-c("RTSE","RNSE","RTSEp")
  rownames(final_table)<-c("LV_Indep","LV_Dep","Hist","Hist_unwt","2011_Avg","Naive","NEST_mcv",
                           "NEST_orc")
 
  save.image("baseball_All_3d.RData")
