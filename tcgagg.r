
#TCGA官网https://portal.gdc.cancer.gov/下载数据，文件夹里准备解压后的Counts文件夹和metadata.json
#数据路径"

library(rjson)
json <- jsonlite::fromJSON("luad-2023-8-16.json")   #metadata?ļ???"
View(json)

sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  

count_file <- list.files('luda-2023-8-16//',pattern = '.tsv',recursive = TRUE)  #Counts?ļ?????"
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})

#???????޸Ļ?????
matrix = data.frame(matrix(nrow=60660,ncol=0))

#???????޸?????????
for (i in 1:557){
    path = paste0('luda-2023-8-16//',count_file[i])   #Counts?ļ?????
    data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
    colnames(data)<-data[2,]
    data <-data[-c(1:6),]
  data <- data[8]   #数据类型，选择其中之一 3：unstranded；4：stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
    matrix <- cbind(matrix,data)
}



 write.csv(matrix,'luad.csv',row.names = TRUE)
