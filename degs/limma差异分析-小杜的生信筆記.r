## limma 做差异分析
## Author: 小杜的生信筆記

setwd("D:\\小杜的生信筆記")

## 导入包
library(limma)
library(dplyr)
## 加载数据
load(file = "luda-limma-8-22-75mut.rda")
df <- log2(df + 1) 
df =  scale(df)

head(df)
### 样本信息注释
list <- c(rep("Wild", 429), rep("Mut",75)) %>% factor(., levels = c("Wild", "Mut"), ordered = F)
head(list)
list <- model.matrix(~factor(list)+0) 
colnames(list) <- c( "Wild","Mut")
df.fit <- lmFit(df, list)

##  差异分析
df.matrix <- makeContrasts(Wild - Mut, levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
head(tempOutput)

## 
nrDEG = na.omit(tempOutput) ## 去掉数据中有NA的行或列
diffsig <- nrDEG  
write.csv(diffsig, "luad-Mut-v-Wild.csv")  


##  筛选出差异基因
foldChange = 1
padj = 0.05
All_diffSig <- diffsig[(diffsig$P.Value < padj & (diffsig$logFC>foldChange | diffsig$logFC < (-foldChange))),]
dim(All_diffSig)
write.csv(All_diffSig, "luad-degs-8-20.csv")  ##输出差异基因数据集


## 筛选上调和下调的差异
diffup <-  All_diffSig[(All_diffSig$adj.P.Val < padj & (All_diffSig$logFC > foldChange)),]
dim(diffup)
write.csv(diffup, "luad-diffup-8-20.csv")
#
diffdown <- All_diffSig[(All_diffSig$adj.P.Val < padj & (All_diffSig$logFC < -foldChange)),]
dim(diffdown)
write.csv(diffdown, "luad-diffdown-8-20.csv")

#---------------------------------
# 绘制火山图
library(ggplot2)
library(ggrepel)
###
logFC <- diffsig$logFC
deg.padj <- diffsig$P.Value
data <- data.frame(logFC = logFC, padj = deg.padj)
data$group[(data$padj > 0.05 | data$padj == "NA") | (data$logFC < foldChange) & data$logFC > -foldChange] <- "Not"
data$group[(data$padj <= 0.05 & data$logFC > 1)] <-  "Up"
data$group[(data$padj <= 0.05 & data$logFC < -1)] <- "Down"
x_lim <- max(logFC,-logFC)
###
pdf('volcano.pdf',width = 7,height = 6.5)
label = subset(diffsig,P.Value <0.05 & abs(logFC) > 1)
label1 = rownames(label)

colnames(diffsig)[1] = 'log2FC'
Significant=ifelse(
  (diffsig$adj.P.Val < 0.05 & abs(diffsig$log2FC)> 1), 
  ifelse(diffsig$log2FC > 1,"Up","Down"), "Not")

p =ggplot(diffsig, aes(log2FC, -log10(P.Value)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#0072B6","grey","#BC3C28"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-1,1), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(p.value)")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  str(diffsig, max.level = c(-1, 1))+theme_bw()

ggsave("火山图-degs-8-20.png",p,width = 5800,height = 3000,units = "px")



##  绘制差异基因热图
library(pheatmap)

## 本章节，我们是只使用count值进行热图绘制，后续的分子中，我们会对其进行优化
DEG_id <- read.csv("luad-degs-8-22.csv", header = T)
head(DEG_id)
## 匹配
DEG_id <- unique(DEG_id$X)
DEG_exp <- df[DEG_id,]
hmexp <- na.omit(DEG_exp)

## 样本注释信息 
annotation_col <- data.frame(Group = factor(c(rep("Wild", 429), rep("Mut",75))))
rownames(annotation_col) <- colnames(hmexp)

##  绘制热图 

p1 <- pheatmap(hmexp,
              annotation_col = annotation_col,
              color = colorRampPalette(c("blue","white","red"))(50),
              cluster_cols = F,
              show_rownames = F,
              show_colnames = F,
              scale = "row", ## none, row, column
              fontsize = 12,
              fontsize_row = 12,
              fontsize_col = 6,
              border = FALSE)
ggsave("热图-degs-8-22.pdf",p1,width = 6,height = 3)
