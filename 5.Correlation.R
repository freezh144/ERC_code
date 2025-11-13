rm(list = ls())
setwd(getwd())
##load library
library(psych)
library(reshape2)
#risk.score
risk<-read.csv("risk_score.csv",header=TRUE,row.names = 1)
head(risk)
dim(risk)

#exp 
library(data.table)
M<-fread("exp.csv",header=TRUE)
M<-as.data.frame(M)
row.names(M)<-M$V1
M<-M[,-1]
head(M)[1:6]

#correlation
cor = corr.test(M1, risk, method="spearman", adjust="fdr")
cmt <-cor$r
head(cmt)
pmt <- cor$p
head(pmt)

cmt<-as.data.frame(cmt)
cmt1<-subset(cmt,cmt$Risk_score>0.3 |  cmt$Risk_score<(-0.3))
dim(cmt1)
head(cmt1)

Correlation<-cmt1
Correlation$Correlation<-ifelse(cmt1$Risk_score>0.3,"Positive_correlation","Negative_correlation" )
Correlation<-Correlation[,-1]
head(Correlation)
names(Correlation)<-row.names(cmt1)

genes<-row.names(cmt1)
head(genes)
M4<-M1[,genes]
head(M4)[,1:6]
dim(M4)


#Enrichment analysis
rm(list = ls())
exp_genes<-read.csv("exp_selcect.csv",header=TRUE,row.names = 1)
head(exp_genes)[,1:6]
genes<-colnames(exp_genes)
head(genes)
genes<-as.data.frame(genes)
colnames(genes)[1]<-"symbol"
head(genes)
library(dplyr)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
#BiocManager::install( "topGO" )
library(topGO)
library(enrichplot)
s2e <- bitr(unique(genes$symbol), fromType = "SYMBOL",  #ID转换核心函数bitr
            toType = c( "ENTREZID"),
            OrgDb = org.Hs.eg.db)
head(s2e)
head(genes)
deg <- inner_join(genes,s2e,by=c("symbol"="SYMBOL"))
head(deg)
dim(deg)
gene_ENTREZID<-deg$ENTREZID

library(ggplot2)
library(ggpubr)
#GO
ego<-enrichGO(
  gene = gene_ENTREZID,   
  keyType = "ENTREZID",  
  OrgDb = org.Hs.eg.db,   
  ont = "all",            
  pAdjustMethod = "BH",   
  pvalueCutoff = 0.05,    
  readable = TRUE
)
head(ego)
write.csv(ego,"relative_genes_ego_output.csv")

##GO dotplot
pdf(file = "GO_dotplot.pdf",width = 8,height = 8)
dotplot(ego,showCategory = 10,orderBy = "x")
dev.off()

##KEGG
ekegg<-enrichKEGG(gene = gene_ENTREZID,
                  keyType = "kegg", 
                  organism  = "hsa",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH")
head(ekegg)
write.csv(ekegg,"dif_genes_ekegg_output.csv")

#KEGG dotplot
pdf(file = "KEGG_dotplot.pdf",width = 9,height = 8)
dotplot(ekegg,showCategory= 10,orderBy = "x")+
  scale_fill_continuous(low="#e06663",high="#327eba",name="p.adjust")
dev.off()

