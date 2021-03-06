
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


#1. Fig 1
load('fig1.RData')


plotdata1<- as.data.frame(c(colMeans(out_m10$risk.gl)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.km)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbsg)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbm)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.ds_mcv)/colMeans(out_m10$risk.ds_orc)))
names(plotdata1)<-"Relative_Risk"
plotdata1$type<-as.factor(rep(c('Grp Linear','NPMLE','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata1$u<- rep(sigma2_u,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=u,y=Relative_Risk,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=u,y=Relative_Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(0.975,1.32))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata2<- as.data.frame(c(colMeans(out_m30$risk.gl)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.km)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbsg)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbm)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.ds_mcv)/colMeans(out_m30$risk.ds_orc)))
names(plotdata2)<-"Risk"
plotdata2$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata2$u<- rep(sigma2_u,5)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(0.975,1.32))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata3<- as.data.frame(c(colMeans(out_m50$risk.gl)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.km)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbsg)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbm)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.ds_mcv)/colMeans(out_m50$risk.ds_orc)))
names(plotdata3)<-"Risk"
plotdata3$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata3$u<- rep(sigma2_u,5)
g3<-ggplot()+geom_line(data=plotdata3,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata3,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(0.975,1.32))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

ggarrange(g1,g2,g3,ncol=3,nrow=1,common.legend = TRUE)
#----------------------------------------------------------
#2. Fig 2
load('fig2.RData')


plotdata1<- as.data.frame(c(colMeans(out_m10$risk.gl)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.km)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbsg)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbm)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.ds_mcv)/colMeans(out_m10$risk.ds_orc)))
names(plotdata1)<-"Relative_Risk"
plotdata1$type<-as.factor(rep(c('Grp Linear','NPMLE','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata1$u<- rep(sigma2_u,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=u,y=Relative_Risk,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=u,y=Relative_Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(1,1.16))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata2<- as.data.frame(c(colMeans(out_m30$risk.gl)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.km)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbsg)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbm)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.ds_mcv)/colMeans(out_m30$risk.ds_orc)))
names(plotdata2)<-"Risk"
plotdata2$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata2$u<- rep(sigma2_u,5)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(1,1.16))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata3<- as.data.frame(c(colMeans(out_m50$risk.gl)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.km)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbsg)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbm)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.ds_mcv)/colMeans(out_m50$risk.ds_orc)))
names(plotdata3)<-"Risk"
plotdata3$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata3$u<- rep(sigma2_u,5)
g3<-ggplot()+geom_line(data=plotdata3,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata3,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(1,1.16))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

ggarrange(g1,g2,g3,ncol=3,nrow=1,common.legend = TRUE)
#--------------------------------------------------------------
#3. Fig 3
load('fig3.RData')


plotdata1<- as.data.frame(c(colMeans(out_m10$risk.gl)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.km)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbsg)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbm)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.ds_mcv)/colMeans(out_m10$risk.ds_orc)))
names(plotdata1)<-"Relative_Risk"
plotdata1$type<-as.factor(rep(c('Grp Linear','NPMLE','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata1$u<- rep(sigma2_u,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=u,y=Relative_Risk,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=u,y=Relative_Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(0.975,1.85))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata2<- as.data.frame(c(colMeans(out_m30$risk.gl)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.km)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbsg)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbm)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.ds_mcv)/colMeans(out_m30$risk.ds_orc)))
names(plotdata2)<-"Risk"
plotdata2$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata2$u<- rep(sigma2_u,5)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(0.975,1.85))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata3<- as.data.frame(c(colMeans(out_m50$risk.gl)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.km)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbsg)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbm)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.ds_mcv)/colMeans(out_m50$risk.ds_orc)))
names(plotdata3)<-"Risk"
plotdata3$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata3$u<- rep(sigma2_u,5)
g3<-ggplot()+geom_line(data=plotdata3,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata3,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(0.975,1.85))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

ggarrange(g1,g2,g3,ncol=3,nrow=1,common.legend = TRUE)
#--------------------------------------------------------------
#4. Fig 4
load('fig4.RData')


plotdata1<- as.data.frame(c(colMeans(out_m10$risk.gl)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.km)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbsg)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbm)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.ds_mcv)/colMeans(out_m10$risk.ds_orc)))
names(plotdata1)<-"Relative_Risk"
plotdata1$type<-as.factor(rep(c('Grp Linear','NPMLE','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata1$u<- rep(sigma2_u,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=u,y=Relative_Risk,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=u,y=Relative_Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(0.975,1.3))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata2<- as.data.frame(c(colMeans(out_m30$risk.gl)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.km)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbsg)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbm)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.ds_mcv)/colMeans(out_m30$risk.ds_orc)))
names(plotdata2)<-"Risk"
plotdata2$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata2$u<- rep(sigma2_u,5)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(0.975,1.3))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata3<- as.data.frame(c(colMeans(out_m50$risk.gl)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.km)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbsg)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbm)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.ds_mcv)/colMeans(out_m50$risk.ds_orc)))
names(plotdata3)<-"Risk"
plotdata3$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata3$u<- rep(sigma2_u,5)
g3<-ggplot()+geom_line(data=plotdata3,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata3,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(0.975,1.3))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

ggarrange(g1,g2,g3,ncol=3,nrow=1,common.legend = TRUE)
#--------------------------------------------------------------
#5. Fig 5
load('fig5.RData')


plotdata1<- as.data.frame(c(colMeans(out_m10$risk.gl)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.km)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbsg)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbm)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.ds_mcv)/colMeans(out_m10$risk.ds_orc)))
names(plotdata1)<-"Relative_Risk"
plotdata1$type<-as.factor(rep(c('Grp Linear','NPMLE','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata1$u<- rep(sigma2_u,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=u,y=Relative_Risk,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=u,y=Relative_Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(1,1.175))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata2<- as.data.frame(c(colMeans(out_m30$risk.gl)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.km)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbsg)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbm)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.ds_mcv)/colMeans(out_m30$risk.ds_orc)))
names(plotdata2)<-"Risk"
plotdata2$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata2$u<- rep(sigma2_u,5)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(1,1.175))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata3<- as.data.frame(c(colMeans(out_m50$risk.gl)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.km)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbsg)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbm)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.ds_mcv)/colMeans(out_m50$risk.ds_orc)))
names(plotdata3)<-"Risk"
plotdata3$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata3$u<- rep(sigma2_u,5)
g3<-ggplot()+geom_line(data=plotdata3,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata3,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(1,1.175))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

ggarrange(g1,g2,g3,ncol=3,nrow=1,common.legend = TRUE)
#--------------------------------------------------------------
#6. Fig 6
load('fig6.RData')


plotdata1<- as.data.frame(c(colMeans(out_m10$risk.gl)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.km)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbsg)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.xkbm)/colMeans(out_m10$risk.ds_orc),
                            colMeans(out_m10$risk.ds_mcv)/colMeans(out_m10$risk.ds_orc)))
names(plotdata1)<-"Relative_Risk"
plotdata1$type<-as.factor(rep(c('Grp Linear','NPMLE','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata1$u<- rep(sigma2_u,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=u,y=Relative_Risk,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=u,y=Relative_Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(1,1.7))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata2<- as.data.frame(c(colMeans(out_m30$risk.gl)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.km)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbsg)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.xkbm)/colMeans(out_m30$risk.ds_orc),
                            colMeans(out_m30$risk.ds_mcv)/colMeans(out_m30$risk.ds_orc)))
names(plotdata2)<-"Risk"
plotdata2$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata2$u<- rep(sigma2_u,5)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(1,1.7))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

plotdata3<- as.data.frame(c(colMeans(out_m50$risk.gl)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.km)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbsg)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.xkbm)/colMeans(out_m50$risk.ds_orc),
                            colMeans(out_m50$risk.ds_mcv)/colMeans(out_m50$risk.ds_orc)))
names(plotdata3)<-"Risk"
plotdata3$type<-as.factor(rep(c('Grp Linear','KM','XKB.SG','XKB.M',
                                'DS'),each=length(sigma2_u)))
plotdata3$u<- rep(sigma2_u,5)
g3<-ggplot()+geom_line(data=plotdata3,aes(x=u,y=Risk,color=type),size=1)+
  geom_point(data=plotdata3,aes(x=u,y=Risk,color=type,shape=type),size=4,fill=NA)+
  theme_bw()+ylab("Relative Risk")+scale_y_continuous(limits = c(1,1.7))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

ggarrange(g1,g2,g3,ncol=3,nrow=1,common.legend = TRUE)



