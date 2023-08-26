#b站 医学生吉克

rm(list = ls())
setwd("D:/R_do/bili")   #设置工作路径
library(tidyverse)
library(GSEABase)
library(GSVA)

load("luda-limma-8-22-75mut.rda")

load("cellMarker_ssGSEA.Rdata")
df <- log2(df + 1) 
expr <- as.matrix(df)   #将expr转换为矩阵格式

#免疫浸润
gsva_data <- gsva(expr,
                  cellMarker, 
                  method = "ssgsea")

gsva_data1 <- gsva_data %>% t() %>% as.data.frame()
write.csv(gsva_data1,'luad-gsva-end.csv',row.names = TRUE)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
cooo <- cor(gsva_data1)

corrplot(cooo, 
         col=colorRampPalette(c('navy', 'white', 'red'))(20),
         tl.cex = 0.6,
         tl.col="black") 
write.csv(cooo,"corr.csv",row.names = TRUE)

blue_palette <- colorRampPalette(c("navy", "white", "red"))(20)
p = pheatmap(cooo,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         cellwidth = 16,
         cellheight = 14,  
         color = blue_palette,
         show_cellnote = TRUE,
         digits = 2,
         angle_col = 90,
         ylab = "right",
         display_numbers = TRUE,
         fontsize = 8,
         colnames = TRUE,
         number_color = "black",
         border=FALSE
         )
ggsave("ssgsea-luad.png",p,width = 4500,height = 2000,units = "px")

#自定义
geneset <- getGmt("Imm_terms.gmt")

es.max <- gsva(expr, 
               geneset, #刚自定义的通路相关基因集
               method='ssgsea', 
               kcdf='Gaussian', 
               abs.ranking=TRUE)










