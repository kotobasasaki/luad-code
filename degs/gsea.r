library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(ggridges)
library(msigdbr)
rt = read.csv(file = "luad-degs-8-22.csv", 
              sep = ",",
              header = T)
genes=as.vector(rt[,1])

#取rt的第一列，即基因名字，将其转换为向量，并赋值给genes变量

gene_map <- select(org.Hs.eg.db, keys=genes, keytype="SYMBOL", columns=c("ENTREZID"))
colnames(gene_map)[1]<-"Gene"
aaa<-inner_join(gene_map,rt,by = "Gene")
aaa<-aaa[,-1]
aaa<-na.omit(aaa)
aaa$logFC<-sort(aaa$logFC,decreasing = T)
geneList = aaa[,2]
names(geneList) = as.character(aaa[,1])
geneList

#GSEA分析——GO
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
head(Go_gseresult)
#GSEA分析——KEGG
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
head(KEGG_gseresult)
#GSEA分析——Reactome
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
head(Go_Reactomeresult)

ridgeplot(Go_gseresult,10)
#富集曲线图类型1：
gseaplot(Go_gseresult,1,pvalue_table = TRUE) #输出第1个结果

#富集曲线图类型2：
gseaplot2(Go_gseresult,100,pvalue_table = TRUE)#输出第212个结果
#gseaplot2还可以同时显示复数个功能组的富集曲线，并标记P值：
gseaplot2(Go_Reactomeresult, 9:12, pvalue_table = TRUE)
