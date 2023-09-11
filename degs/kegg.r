library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)
R.utils::setOption("clusterProfiler.download.method",'auto')

rt = read.csv(file = "luad-degs-8-22.csv", 
              sep = ",",
              header = T)
genes=as.vector(rt[,1])

#取rt的第一列，即基因名字，将其转换为向量，并赋值给genes变量
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)

out=cbind(rt,entrezID=entrezIDs)
gene=out$entrezID  
kk2<- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05) 
KEGG = data.frame(kk2)


write.csv(KEGG,'luad-kegg-degs-8-22.csv')

kegg1=data.frame(Category = "ALL",
                 ID = KEGG$ID,
                 Term = KEGG$Description, 
                 Genes = gsub("/", ", ", KEGG$geneID), 
                 adj_pval = KEGG$p.adjust)

p = dotplot(kk2, showCategory = 3)
ggsave("kegg-degs-8-22.png",p,width = 5000,height = 2000,units = "px")
genelist <- data.frame(ID = out$entrezID, logFC = rt$logFC)

row.names(genelist)=genelist[,1]
circ <- circle_dat(kegg1, genelist)


GOBubble(circ, labels = 3,table.legend =F)

GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=3) 

termNum =3 #限定term数目

geneNum = nrow(genelist) #限定基因数目

chord2 <- chord_dat(circ, genelist[1:geneNum,], kegg1$Term[1:termNum])



p1 = GOHeat(chord2, nlfc =1, fill.col = c('red', 'white', 'blue'))

ggsave("kegg-degs-heat-8-22.png",p1,width = 5000,height = 2000,units = "px")
