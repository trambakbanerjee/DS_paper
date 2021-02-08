

require(ggplot2)
require(gridExtra)
library(ggpubr)


load('fig11_left_mi.RData')
g4<-g1
load('fig11_center_mi.RData')
g5<-g1
load('fig11_right_mi.RData')
g7<-g1
load('fig12_left_mi.RData')
g9<-g1
load('fig12_center_mi.RData')
g12<-g1
load('fig12_right_mi.RData')
g13<-g1

ggarrange(g4,g5,g7,ncol=3,nrow=1,common.legend = TRUE)
ggarrange(g9,g12,g13,ncol=3,nrow=1,common.legend = TRUE)


