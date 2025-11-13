rm(list = ls())
setwd(getwd())
library(survival)   
library(survminer)
##input data
expr_cli<-fread("data.csv",header = T)
expr_cli<-as.data.frame(expr_cli)
row.names(expr_cli)<-expr_cli$V1
expr_cli<-expr_cli[,-1]
head(expr_cli)[,1:10]
multi_expr_cli<-expr_cli
head(multi_expr_cli)[,1:10]
dim(multi_expr_cli)
covariates<-colnames(multi_expr_cli)[3:5]
multi_expr_cli<-multi_expr_cli[,c("Status","Time",covariates   )]
head(multi_expr_cli)
res.cut <- surv_cutpoint(multi_expr_cli, time = "Time", event = "Status",
                         variables = covariates )
res.cut

for (i in 1:ncol(multi_expr_cli)) {
  n=colnames(multi_expr_cli)[i]
  multi_expr_cli[,n]=ifelse(multi_expr_cli[,n]>=res.cut[n,1],"High","Low")
}
head(multi_expr_cli)

multi_expr_cli_high_low<-cbind(expr_cli[,1:2],multi_expr_cli[,3:ncol(multi_expr_cli)])
head(multi_expr_cli_high_low)

#KM picutres
genes<-names(multi_expr_cli_high_low)[3:ncol(multi_expr_cli_high_low)]
head(genes)
length(genes)
plist<-list()
for (i in genes) {
  print(i)
  fit<-survfit(Surv(Time,Status)~multi_expr_cli_high_low[,i],data=multi_expr_cli_high_low)
  p<-ggsurvplot(fit,linetype = "strata",
                pval = TRUE,
                title=i,
                palette = c("red","blue" ),legend = c(0.8, 0.2),
                legend.labs=c(paste("High"),paste("Low"))
  )
  png(paste(i,"_surv.png"))
  dev.off()
}

