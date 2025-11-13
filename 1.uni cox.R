#Single gene Cox analysis
rm(list = ls())
setwd(getwd())
#input data
library(data.table)
expr_cli<-fread("Exp.csv",header = T)
expr_cli<-as.data.frame(expr_cli)
row.names(expr_cli)<-expr_cli$V1
expr_cli<-expr_cli[,-1]
head(expr_cli)[,1:6]

#cox
#load library
library(survival)  
library(survminer)
covariates<-names(expr_cli)[3:ncol(expr_cli)]  
univ_formulas<-sapply(covariates,  function(x) 
  as.formula(paste("Surv(os_followup,os_status)~",x)) )     
head(univ_formulas)
univ_models<-lapply(univ_formulas,function(x){coxph(x,data = expr_cli)})   
head(univ_models)

#HR,95%CI and p-value
univ_results<-lapply(univ_models, function(x) {
  x<-summary(x)
  p.value<-signif(x$wald["pvalue"],digits = 2) 
  HR<-signif(x$coef[2],digits = 2)   
  HR.confint.lower<-signif(x$conf.int[,"lower .95"],2)   
  HR.confint.upper<-signif(x$conf.int[,"upper .95"],2) 
  HR<-paste(HR,"(",HR.confint.lower,"-",HR.confint.upper,")"    )
  res<-c(p.value,HR)
  names(res)<-c("p.value","HR (95% CI for HR)")
  return(res)
})
head(univ_results)
res<-t(as.data.frame(univ_results,check.names=FALSE))
res<-as.data.frame(res)
head(res)

res_sig1<-res[which(res$p.value <0.05), ]
head(res_sig1)
dim(res_sig1)
write.csv(res_sig1,"df_univariate_cox_result.csv")



