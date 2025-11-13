rm(list = ls())
setwd(getwd())
#load library
library(glmnet)
library(survivalROC)
library(plotROC)
library(survival)   
library(survminer)

#input data
res_cox<-read.csv("multi_expr.csv",header = T,row.names = 1)
head(res_cox)#[,1:6]
dim(res_cox)

#cox
multi_res.cox<-coxph(Surv(os_followup,os_status)~.,data = res_cox)
x<-summary(multi_res.cox)
coefficients<-multi_res.cox$coefficients  
coefficients
length(coefficients)
pvalue<-signif(as.matrix(x$coefficients)[,5],2)
head(pvalue)
HR<-signif(as.matrix(x$coefficients)[,2],2)
head(HR)
low<-signif(x$conf.int[,3],2)
high<-signif(x$conf.int[,4],2)
multi_res<-data.frame(p.value=pvalue,
                      HR=paste(HR,"(",low,"-",high,")",sep = ""),
                      stringsAsFactors = F)

##forest plot
pdf("Hazard ratios of multi genes.pdf",width = 10,height = 10)
ggforest(model = multi_res.cox,data = multi_expr_cli,
         main = "Hazard ratios of multi genes",fontsize = 1)
dev.off()




