rm(list = ls()) 
setwd(getwd())
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

library(data.table)
testExpr<-fread("Exp.csv",header = T )
testExpr<-as.data.frame(testExpr)
head(testExpr)[,1:6]

#Calculate
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )


#Plot
drug<- read.csv("drug_select.csv",header=TRUE,row.names = 1)
head(drug)[,1:5]
p<-ggplot(drug,aes(x=Risk_score, y=Daporinad,fill=Risk_score,
                   color=Risk_score))+
  geom_rain(alpha=0.6,size=0.5)+
  theme_pubr(legend = "bottom")+ ggsci::scale_fill_d3(alpha = 0.6)+
  ggsci::scale_color_d3(alpha = 0.6)+
  ylab("Estimated IC50")+
  stat_compare_means(aes(x=Risk_score,y=Daporinad,group=Risk_score),
                     method="t.test",show.legend = F,
                     label = "p.format",  ,hide.ns = TRUE)
p




