rm(list = ls())
setwd(getwd())
#load library
library(IOBR)
library(EPIC)
library(data.table)
eset_stad<-fread( "data.csv",header = T)
eset_stad<-as.data.frame(eset_stad)
eset_stad[1:5, 1:5]

#1. TIMER
timer<-deconvo_tme(eset = eset_stad, method = "timer", 
                   group_list = rep("stad",dim(eset_stad)[2]))
head(timer)
res<-cell_bar_plot(input = timer, features = colnames(timer),
                    title = "TIMER Cell Fraction")
res


#2. Boxplot
dat1<-read.csv("dat1.csv",header=TRUE   )
head(dat1)
col_pal<-c("red","#04DBB4" )
p1<-ggplot(dat1,aes(Cell_type,Proportion,fill = Risk_score)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_pubr() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "bottom") + coord_flip()+
  scale_fill_manual(values = col_pal)+ 
  stat_compare_means(aes(group = Risk_score,label = ..p.signif..),
                     method = "kruskal.test")
print(P1)

##3. Correlation
gene<-read.csv("multi_expr.csv",header=TRUE,row.names = 1   )
head(gene)
head(cell)[,1:6]
library(linkET)
library(vegan)
mantel <- mantel_test(gene,cell, 
                      spec_select = list(IL20RB = 1,
                                         NR1H2=2,
                                         SSTR5=3
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), # 
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), 
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
p2<-qcorrplot(correlate(cell), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd),
              data = mantel, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
print(P2)


