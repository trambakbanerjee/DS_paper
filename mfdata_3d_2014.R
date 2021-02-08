
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
library(readr)

#1. Get monthly treasury yields from 2014
treasury_3month <- read_csv("treasury_3month.csv", 
                            col_types = cols(date = col_date(format = "%m/%d/%Y")))

treasury_2015<- treasury_3month[treasury_3month$date>="2014-01-01",]
date.unique<-unique(treasury_2015$date)
tr_2015<-as.data.frame(date.unique)
names(tr_2015)<-"date"
tr_2015$threemnth <- matrix(0,length(date.unique))
tr_2015$ret<- matrix(0,length(date.unique))
for(i in 1:length(date.unique)){
  
  tr_2015$threemnth[i]<- mean(treasury_2015$`3mo_treasury`[treasury_2015$date==tr_2015$date[i]])
  tr_2015$ret[i]<- tr_2015$threemnth[i]/100
 }
tr_2015$ym<- format(tr_2015$date,"%Y-%m")

mfdata <- read_csv("finance data/mutualfunds/mfdata.csv", 
                   col_types = cols(caldt = col_date(format = "%Y%m%d"), 
                  crsp_fundno = col_integer()))
mfdata<-mfdata[complete.cases(mfdata$mnav),]
mfdata<-mfdata[complete.cases(mfdata$mret),]
mfdata<-mfdata[complete.cases(mfdata$mtna),]
mfdata<-mfdata[mfdata$mtna>0,]
mfdata<-mfdata[complete.cases(mfdata$crsp_fundno),]

mfdata_part1<-mfdata[mfdata$caldt>="2014-01-31" & mfdata$caldt<="2014-06-30",]
mfdata_part1$ym<- format(mfdata_part1$caldt,"%Y-%m")
cc<-as.data.frame(table(mfdata_part1$crsp_fundno))
dd<- as.data.frame(cc$Var1[cc$Freq>=4])
names(dd)<-'crsp_fundno'
idx<-which(mfdata_part1$crsp_fundno %in% dd$crsp_fundno)
mfdata_part1<-as.data.frame(mfdata_part1[idx,])
temp<-as.data.frame(merge(x = mfdata_part1, y = tr_2015, by = "ym", all.x = TRUE))
idx<-is.na(temp$crsp_fundno)
temp<-temp[!idx,]
mfdata_part1<- temp[,c("crsp_fundno","ym","mret","ret")]
n<-length(unique(mfdata_part1$crsp_fundno))
funds_part1<- unique(mfdata_part1$crsp_fundno)
set.seed(1)
idx<-sample(1:n,5000)
funds_part1<-funds_part1[idx]
temp<-mfdata_part1[mfdata_part1$crsp_fundno %in% funds_part1,]
mfdata_part1<-temp

mfdata_part2<-mfdata[mfdata$caldt>="2014-07-31" & mfdata$caldt<="2014-12-31",]
idd<-which(mfdata_part2$crsp_fundno %in% mfdata_part1$crsp_fundno)
mfdata_part2<-as.data.frame(mfdata_part2[idd,])
cc<-as.data.frame(table(mfdata_part2$crsp_fundno))
dd<- as.data.frame(cc$Var1[cc$Freq>=3])
names(dd)<-'crsp_fundno'
idx<-which(mfdata_part2$crsp_fundno %in% dd$crsp_fundno)
mfdata_part2<-as.data.frame(mfdata_part2[idx,])
mfdata_part2$ym<- format(mfdata_part2$caldt,"%Y-%m")
temp<-as.data.frame(merge(x = mfdata_part2, y = tr_2015, by = "ym", all.x = TRUE))
idx<-is.na(temp$crsp_fundno)
temp<-temp[!idx,]
mfdata_part2<- temp[,c("crsp_fundno","ym","mret","ret")]


#----------------------------------
n<-length(unique(mfdata_part1$crsp_fundno))
funds_part1<- unique(mfdata_part1$crsp_fundno)

X<-list()
S2.sample<-matrix(0,n,1)
m<-matrix(0,n,1)
for(i in 1:n){
  
  idx<-which(mfdata_part1$crsp_fundno==funds_part1[i])
  m[i]<-length(idx)
  X[[i]]<- (mfdata_part1$mret[idx]-mfdata_part1$ret[idx])
  S2.sample[i]<-var(X[[i]])
}
S.sample<-sqrt(S2.sample)
xbar<-sapply(1:n,function(i) mean(X[[i]]))

funds_part2<- unique(mfdata_part2$crsp_fundno)
idy<- sapply(1:length(funds_part2),
             function(i) which(funds_part2[i] == funds_part1))

theta<-matrix(0,length(funds_part2),1)
mm<-matrix(0,length(funds_part2),1)
for(i in 1:length(funds_part2)){
  
  idx<-which(mfdata_part2$crsp_fundno==funds_part2[i])
  mm[i]<- length(idx)
  temp<-(mfdata_part2$mret[idx]-mfdata_part2$ret[idx])
  theta[i]<- mean(temp)/sd(temp)
}
rm('temp')

#--------------------- Methods ------------------------
nn<-length(theta)
#1. Group Linear
thetahat.gl<- grouplinear(xbar,S2.sample/m)
risk2.gl_r<- sum((theta-(thetahat.gl[idy]/S.sample[idy]))^2)/nn
risk1.gl_r<- sum(mm*(theta-(thetahat.gl[idy]/S.sample[idy]))^2)/nn
risk0.gl_r<- sum(mm*abs(1-(thetahat.gl[idy]/S.sample[idy])/theta))/nn

#2. XKB - 3 versions
thetahat.xkb.g<- thetahat.G(xbar,S2.sample/m)
risk2.xkb_g<- sum((theta-(thetahat.xkb.g[idy]/S.sample[idy]))^2)/nn
risk1.xkb_g<- sum(mm*(theta-(thetahat.xkb.g[idy]/S.sample[idy]))^2)/nn
risk0.xkb_g<- sum(mm*abs(1-(thetahat.xkb.g[idy]/S.sample[idy])/theta))/nn

thetahat.xkb.sg<- thetahat.SG(xbar,S2.sample/m)
risk2.xkb_sg<- sum((theta-(thetahat.xkb.sg[idy]/S.sample[idy]))^2)/nn
risk1.xkb_sg<- sum(mm*(theta-(thetahat.xkb.sg[idy]/S.sample[idy]))^2)/nn
risk0.xkb_sg<- sum(mm*abs(1-(thetahat.xkb.sg[idy]/S.sample[idy])/theta))/nn

thetahat.xkb.m<- thetahat.M(xbar,S2.sample/m)
risk2.xkb_m<- sum((theta-(thetahat.xkb.m[idy]/S.sample[idy]))^2)/nn
risk1.xkb_m<- sum(mm*(theta-(thetahat.xkb.m[idy]/S.sample[idy]))^2)/nn
risk0.xkb_m<- sum(mm*abs(1-(thetahat.xkb.m[idy]/S.sample[idy])/theta))/nn

#3. KM #Does not converge#
X.vec<-unlist(X)
id.km<-rep(1:n,m)
#km<- tryCatch(WGLVmix(X.vec,id.km,rep(1,length(X.vec)),rtol=1e-10),
#              error=function(e) NULL)
# km$dy[km$dy<=0]<-S2.sample[km$dy<=0]
# thetahat.km<- km$dx/sqrt(km$dy)
# risk2.km_dep<- sum((theta-thetahat.km[idy])^2)/nn
# risk1.km_dep<- sum(mm*(theta-thetahat.km[idy])^2)/nn
# risk0.km_dep<- sum(mm*abs(1-thetahat.km[idy]/theta))/nn

#km_ind<- tryCatch(WTLVmix(X.vec,id.km,rep(1,length(X.vec)),rtol=1e-10),
#               error=function(e) NULL)
# km_ind$dy[km_ind$dy<=0]<-S2.sample[km_ind$dy<=0]
# thetahat.km<- km_ind$dx/sqrt(km_ind$dy)
# risk2.km_ind<- sum((theta-thetahat.km[idy])^2)/nn
# risk1.km_ind<- sum(mm*(theta-thetahat.km[idy])^2)/nn
# risk0.km_ind<- sum(mm*abs(1-thetahat.km[idy]/theta))/nn
# gc()

#4. Naive
risk2.naive<- sum((theta-(xbar[idy]/S.sample[idy]))^2)/nn
risk1.naive<- sum(mm*(theta-(xbar[idy]/S.sample[idy]))^2)/nn
risk0.naive<- sum(mm*abs(1-(xbar[idy]/S.sample[idy])/theta))/nn

#7. Double Shrinkage  
eta.grid<- c(0.025,0.05,0.075,0.1,0.15,0.2)

#Sample splitting
alpha<-0.5

eps<-lapply(1:n,function(i) norm.rand(m[i],0,S.sample[i],10*i))
ubar<-sapply(1:n,function(i) mean(X[[i]]+alpha*eps[[i]]))
ebar<-sapply(1:n,function(i) mean(eps[[i]]))
S2.u<- sapply(1:n,function(i)var(X[[i]]+alpha*eps[[i]]))
vbar<-sapply(1:n,function(i) mean(X[[i]]-(1/alpha)*eps[[i]]))
S2.v<- sapply(1:n,function(i)var(X[[i]]-(1/alpha)*eps[[i]]))

# Mahalanobis distance
dat<- cbind(scale(xbar),scale(S2.sample),scale(m))
dat_u<- cbind(scale(ubar),scale(S2.u),scale(m))
Omegainv<- cov(dat)
Omega<-spdinv(Omegainv)
Omegainv_u<- cov(dat_u)
Omega_u<-spdinv(Omegainv_u)
K_md<-matrix(0,n,n)
K_md_u<-matrix(0,n,n)
gc()
for(i in 1:n){
  d<-mahala(dat,mu = dat[i,],sigma = Omegainv)
  K_md[,i]<- exp(-0.5*d)
  d<-mahala(dat_u,mu = dat_u[i,],sigma = Omegainv_u)
  K_md_u[,i]<- exp(-0.5*d)
  print(i)
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
  tauhat[tauhat<0]<-(m[tauhat<0]-3)/((m[tauhat<0]-1)*S2.sample[tauhat<0])
  deltaDS<- xbar+(1/(m*tauhat))*rho1_2[,jj]
  loss_true_md[jj]<-sum((theta-(deltaDS[idy]*sqrt(tauhat[idy])))^2)
  
  rho1_4[,jj]<-out_4[,1]
  rho2_4[,jj]<-out_4[,2]
  
  tauhat<- (m-3)/((m-1)*S2.u)-(2/(m-1))*rho2_4[,jj]
  tauhat[tauhat<0]<-(m[tauhat<0]-3)/((m[tauhat<0]-1)*S2.sample[tauhat<0])
  deltaDS<- ubar+(1/(m*tauhat))*rho1_4[,jj]
  loss_mcv_md[jj]<-sum((vbar-deltaDS)^2)
  print(jj)
  
}

ii<- min(which(loss_true_md==min(loss_true_md,na.rm = T)))
eta_orc_md<- eta.grid[ii]

tauhat<-((m-3))/((m-1)*S2.sample)-(2/(m-1))*rho2_2[,ii]
tauhat_orc<- tauhat
deltaDS_orc<- xbar+(1/(m*tauhat_orc))*rho1_2[,ii]
risk2.ds_orc<- sum((theta-(deltaDS_orc[idy]*sqrt(tauhat_orc[idy])))^2)/nn
risk1.ds_orc<- sum(mm*(theta-(deltaDS_orc[idy]*sqrt(tauhat_orc[idy])))^2)/nn
risk0.ds_orc<- sum(mm*abs(1-(deltaDS_orc[idy]*sqrt(tauhat_orc[idy]))/theta))/nn


ii<- min(which(loss_mcv_md==min(loss_mcv_md,na.rm = T)))
eta_mcv_md<- eta.grid[ii]
tauhat<- (m-3)/((m-1)*S2.sample)-(2/(m-1))*rho2_2[,ii]
tauhatDS<- tauhat
deltaDS<- xbar+(1/(m*tauhatDS))*rho1_2[,ii]
risk2.ds<- sum((theta-(deltaDS[idy]*sqrt(tauhatDS[idy])))^2)/nn
risk1.ds<- sum(mm*(theta-(deltaDS[idy]*sqrt(tauhatDS[idy])))^2)/nn
risk0.ds<- sum(mm*abs(1-(deltaDS[idy]*sqrt(tauhatDS[idy]))/theta))/nn

save.image("mfdata_3d_2014data.RData")
