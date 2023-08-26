library(ggplot2)
library(reshape2)

df = read.csv("luad-gsva-end.csv",sep=",",encoding="UTF-8")
df = melt(df)                    # melt出自reshape2包
head(df)                         # 查看转换完成的数据的前几行
p = ggplot(df, aes(x=factor(Mixture,levels =unique(Mixture)),  # 转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
                   y=value, 
                   fill=factor(variable,levels = unique(variable)), 
))+
  labs(
    x="",   # 调整x轴名称
    y="",   # 调整y轴名称
    fill="" # 调整图例名称
  )+geom_bar(position="fill",
stat="identity",colour="black")+theme(axis.text.y = element_text(size = 14))+ 
  theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 13))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 80,vjust = 1,hjust = 1,size = 5,face=2))+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  ylab("Relative Percent")+
  theme(legend.text = element_text(size = 12))

ggsave("gsva-bar-8-22.png",p,width = 8500,height = 4500,units = "px")

#实际上就很简单，把相应的参数放进槽里面即可

