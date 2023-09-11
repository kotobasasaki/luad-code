library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("DOSE")
library("ggnewscale")
library(GOplot)

rt = read.csv(file = "luad-degs-8-22.csv", 
              sep = ",",
              header = T)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)

out=cbind(rt,entrezID=entrezIDs)
gene=out$entrezID 

erich.do<-DOSE::enrichDO(gene=gene,
                         
                         ont = "DO",
                         
                         pvalueCutoff = 0.05,
                         
                         qvalueCutoff = 0.05,
                         
                         readable = T)

do = data.frame(erich.do)

write.csv(do,'luad-do-8-22.csv')


p =barplot(erich.do,showCategory = 30)
ggsave("do-degs-heat-8-20.png",p,width = 5000,height = 4000,units = "px")

p1 = dotplot(erich.do,,showCategory = 25)
ggsave("do-degs-dot-8-22.png",p1,width = 4500,height = 5500,units = "px")

p3 =cnetplot(erich.do,showCategory = 79)
ggsave("do-degs-cnet-8-22.png",p3,width = 5000,height = 5500,units = "px")

do1=data.frame(Category = "ALL",
               ID = do$ID,
               Term = do$Description, 
               Genes = gsub("/", ", ", do$geneID), 
               adj_pval = do$p.adjust)
genelist <- data.frame(ID = rt$X, logFC = rt$logFC)
circ <- circle_dat(do1, genelist)
GOBubble(circ, labels = 3,table.legend =F)
GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10) 

termNum =79 #限定term数目

geneNum = nrow(genelist) #限定基因数目

chord2 <- chord_dat(circ, genelist[1:geneNum,], do1$Term[1:termNum])
p2 = GOHeat(chord2, nlfc =1, fill.col = c('red', 'white', 'blue'))
ggsave("do-degs-heat-8-20.png",p2,width = 4500,height = 6000,units = "px")
