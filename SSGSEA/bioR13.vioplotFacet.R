####公众号:医学学霸帮####

####公众号:医学学霸帮####


#install.packages("ggplot2")

library(ggplot2)
library(reshape2)
      
setwd("D:\\biowolf\\bioR\\13.vioplotFacet")      

#读取输入文件
rt = read.csv(file = "luad-gsva-end.csv", 
              sep = ",",
              header = T,
              check.names=F,
              row.names=1)

x=colnames(rt)[1]
colnames(rt)[1]="Type"


#差异分析
geneSig=c("")
for(gene in colnames(rt)[2:ncol(rt)]){
	rt1=rt[,c(gene,"Type")]
	colnames(rt1)=c("expression","Type")
	p=1
	if(length(levels(factor(rt1$Type)))>2){
		test=kruskal.test(expression ~ Type, data = rt1)
		p=test$p.value
	}else{
		test=wilcox.test(expression ~ Type, data = rt1)
		p=test$p.value
	}
	Sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","")))
	geneSig=c(geneSig,Sig)
}
colnames(rt)=paste0(colnames(rt),geneSig)


#把数据转换成ggplot2文件
data=melt(rt,id.vars=c("Type"))
colnames(data)=c("Type","Gene","Fraction")

genes <- unique(data$Gene)
top_genes <- head(genes, 10)
other_genes <- setdiff(genes, top_genes)




#绘制
p1=ggplot(data, aes(x = Type, y = Fraction, fill = Type)) +
  guides(fill = guide_legend(title = x)) +
  labs(x = x, y = "Fraction") +
  geom_violin() + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) + 
  facet_wrap(~Gene, nrow = 3) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12)) 

ggsave("vio-gsea-8-22.png",p1,width = 6000,height = 3500,units = "px")

#把数据转换成ggplot2文件
#差异分析
#绘制
#输出
