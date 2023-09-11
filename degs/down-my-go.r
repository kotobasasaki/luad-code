library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)


rt = read.csv(file = "luad-degs-8-22.csv", 
                 sep = ",",
                 header = T)
genes=as.vector(rt[,1])

#取rt的第一列，即基因名字，将其转换为向量，并赋值给genes变量
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)

out=cbind(rt,entrezID=entrezIDs)
gene=out$entrezID 
kk <- enrichGO(gene = gene,
               
               OrgDb = org.Hs.eg.db,
               
               pvalueCutoff =0.05,
               
               qvalueCutoff = 0.05,
               
               ont="all",
               
               readable =T)
write.csv(kk, "go-imm-v-degs-8-22.csv") 
go = data.frame(kk)
go1=data.frame(Category = go$ONTOLOGY,
               ID = go$ID,
               Term = go$Description, 
               Genes = gsub("/", ", ", go$geneID), 
               adj_pval = go$p.adjust)


p=dotplot(kk,showCategory = 10,split="ONTOLOGY",font.size=8) + facet_grid(ONTOLOGY~., scale='free')
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')

ggsave("do-degs-8-22.png",p,width = 5500,height = 3500,units = "px")

genelist <- data.frame(ID = rt$X, logFC = rt$logFC)
row.names(genelist)=genelist[,1]
circ <- circle_dat(go1, genelist)

GOBubble(circ, labels = 3,table.legend =F)
GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10)   
termNum = 30 #限定term数目

geneNum = nrow(genelist) #限定基因数目

chord <- chord_dat(circ, genelist[1:geneNum,], go1$Term[1:termNum])



p1 = GOHeat(chord, nlfc =1, fill.col = c('red', 'white', 'blue'))

ggsave("do-degs-heat-8-20.png",p1,width = 5000,height = 5500,units = "px")

