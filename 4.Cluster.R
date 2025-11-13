rm(list=ls())
setwd(getwd())
dat<-read.csv("multi_expr.csv",header = T,row.names = 1)
head(dat)#[,1:6]
colnames(dat)
dat<-t(dat[,-1:-2])
head(dat)[,1:6]

#1. consensus clustering
dat <- sweep(dat,1, apply(dat,1,median,na.rm=T))
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(as.matrix(dat),maxK=6,reps=1000,pItem=0.8,
                               pFeature=1,
                               title="output",
                               clusterAlg="pam",distance="pearson",
                               seed=1234,plot="png",writeTable=F)
results[[3]][["consensusMatrix"]][1:5,1:5]


##2. KM plot
multi_expr_cli_cluster<-read.csv("multi_expr_cli_cluster.csv",header = T,row.names = 1)
head(multi_expr_cli_cluster)
dim(multi_expr_cli_cluster)

library(ggplot2)
library(survminer)
library(survival)   
cluster.cox<-coxph(Surv(os_followup,os_status)~Cluster,data = multi_expr_cli_cluster)
x<-summary(cluster.cox)
x
fit_cluster<-survfit(Surv(os_followup,os_status)~Cluster,data=multi_expr_cli_cluster)
fit_cluster

col_pal<-c("#006e65","#f05fa6" )
p_cluster<-ggsurvplot(fit_cluster,linetype = "strata",
                      pval = TRUE,risk.table = T,conf.int = F,
                      palette = col_pal,xlab="Time (year)" ,
                      legend.labs=c("Cluster_1","Cluster_2")  #,"Cluster_3"
)
print(p_cluster)


##3. boxplot
library(ggplot2)  
library(ggpubr)
library(reshape2)
sig_exp<-read.csv("dat1.csv",header = T,row.names = 1)
head(sig_exp)
mycomprison<-c("Cluster_1","Cluster_2")
col_pal<-c("#006e65","#f05fa6" )
pp<-ggplot(sig_exp,aes(x=gene,y=value,fill=Cluster))+
  geom_boxplot(outlier.size = 1.0,outlier.colour = "white")+ 
  coord_flip()+
  scale_fill_manual(values =col_pal  )+
  xlab("")+ggtitle("")+ylab("Abundance")+
  stat_compare_means(aes(x=gene,y=value,group=Cluster),
                     label = "p.signif",hide.ns = TRUE,
                     method = "wilcox.test")
print(pp)


#4. UMAP
#input data
exp<-read.csv("exp_cluster.csv",row.names = 1,header = TRUE)  
dim(exp)
head(exp)[,1:ncol(exp)]

#install.packages("umap")
library(umap)  
iris.umap = umap::umap(exp[,5:ncol(exp)])
iris.umap
iris.umap$layout

#ruslut
b1<-iris.umap$layout
head(b1)
b1<-as.data.frame(b1)
#b2<-merge(b1,B1,by="row.names")
b2<-cbind(b1,exp[,1:4])
head(b2)
colnames(b2)[1:2]<-c("UMAP_1","UMAP_2")
head(b2)

##plot
library(ggplot2)
library(ggpubr)
#Cluster2
u2<-ggplot(data=b2,aes(x = UMAP_1, y =UMAP_2, color = Cluster_2,shape=Cluster_2)) +
  geom_point(fill="gray",size=2)+theme_pubr()+
  stat_ellipse()+xlab("UMAP_1")+ylab("UMAP_2")+
  ggsci::scale_color_d3(alpha = 0.8)+
  #scale_shape_manual(values = c(15,16,18))+
  labs(title = "K=2")+
  geom_hline(yintercept=0,linetype=4)+geom_vline(xintercept=0,linetype=4)+
  #scale_color_manual(values=colors)+
  theme(legend.position = c(0.9,0.95),
        legend.title=element_blank())   #,legend.position = 'none'
u2

#Cluster3
u3<-ggplot(data=b2,aes(x = UMAP_1, y =UMAP_2, color = Cluster_3,shape=Cluster_3)) +
  geom_point(fill="gray",size=2)+theme_pubr()+
  stat_ellipse()+xlab("UMAP_1")+ylab("UMAP_2")+
  ggsci::scale_color_d3(alpha = 0.8)+
  #scale_shape_manual(values = c(15,16,18))+
  labs(title = "K=3")+
  geom_hline(yintercept=0,linetype=4)+geom_vline(xintercept=0,linetype=4)+
  #scale_color_manual(values=colors)+
  theme(legend.position = c(0.9,0.9),
        legend.title=element_blank())   #,legend.position = 'none'
u3

#Cluster4
u4<-ggplot(data=b2,aes(x = UMAP_1, y =UMAP_2, color = Cluster_4,shape=Cluster_4)) +
  geom_point(fill="gray",size=2)+theme_pubr()+
  stat_ellipse()+xlab("UMAP_1")+ylab("UMAP_2")+
  ggsci::scale_color_d3(alpha = 0.8)+
  #scale_shape_manual(values = c(15,16,18))+
  labs(title = "K=4")+
  geom_hline(yintercept=0,linetype=4)+geom_vline(xintercept=0,linetype=4)+
  #scale_color_manual(values=colors)+
  theme(legend.position = c(0.9,0.9),
        legend.title=element_blank())   #,legend.position = 'none'
u4

#Cluster5
u5<-ggplot(data=b2,aes(x = UMAP_1, y =UMAP_2, color = Cluster_5,shape=Cluster_5)) +
  geom_point(fill="gray",size=2)+theme_pubr()+
  stat_ellipse()+xlab("UMAP_1")+ylab("UMAP_2")+
  ggsci::scale_color_d3(alpha = 0.8)+
  #scale_shape_manual(values = c(15,16,18))+
  labs(title = "K=5")+
  geom_hline(yintercept=0,linetype=4)+geom_vline(xintercept=0,linetype=4)+
  #scale_color_manual(values=colors)+
  theme(legend.position = c(0.9,0.9),
        legend.title=element_blank())   #,legend.position = 'none'
u5

