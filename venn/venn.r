
# 加载R包，没有安装请先安装  install.packages("包名")
library(venn)         #韦恩图（venn 包，适用样本数 2-7）
library(VennDiagram) 


data = read.table("imm-v-degs.txt",sep="\t",header = T,encoding="UTF-8")


venn_list <- list(data[,1], data[,2])  
venn(venn_list,
     zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
     opacity = 0.3,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 0.5,     # 数字大小
     sncs = 1        # 组名字体大小
)
venn.diagram(x=venn_list,
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c('#000000','#000000'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#0000FF','#CD2626'), # 填充色 配色https://www.58pic.com/
             
             category.names = c("DEGs", "Immune") , #标签名
             
             cat.dist = 0.02, # 标签距离圆圈的远近
             
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='两组.png',# 文件保存
             
             imagetype="png",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
             compression = "lzw"# 压缩算法
             
)
venn.diagram(x=venn_list,
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c('#000000','#000000',"#000000"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#FFFF00','#FF0000',"#00BFFF"), # 填充色 配色https://www.58pic.com/
             
             category.names = c("Spe to november 2022", "December 2022 to Feb 2023","Mar to May 2023") , #标签名
             
             cat.dist = 0.02, # 标签距离圆圈的远近
             
             cat.pos = c(-10, 10, -180), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 1.5, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='total.png',# 文件保存
             
             imagetype="png",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
             compression = "lzw"# 压缩算法
             
)

grid.draw(data)
df_inter <- get.venn.partitions(df)
for (i in 1:nrow(df_inter)) df_inter[i,'values'] <- paste(df_inter[[i,'..values..']], collapse = ', ')
df_inter[-c(5, 6)]
write.table(df_inter, 'df_Venn.txt', row.names = FALSE, sep = '\t', quote = FALSE)
