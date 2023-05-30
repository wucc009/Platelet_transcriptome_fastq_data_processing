# ------------------------------------
#       一、导入数据及预处理
# ------------------------------------
setwd("C:\\Users\\13395\\Desktop")
filename = c(list.files(pattern = "txt$"))
data = NULL
for(i in 1:length(filename)) {
  dataname = paste0("data_", formatC(i, width = 4, flag = "0"))
  assign(dataname, read.table(filename[i], header  = T, sep = "\t", comment.char = "")) # assign将读取的数据赋给dataname这个变量下的字符串
  temp = as.matrix(get(dataname)[,c(1,10,9)])
  data = rbind(data, temp) # 取出需要的列按行合并为一个总得矩阵
}
data = as.data.frame(unique(data)) # 去冗余
data$gene_id = substr(data$gene_id, start = 1, stop = 15)
# ------------------------------------
#          二、添加circBankID
# ------------------------------------
setwd("C:\\Users\\13395\\Desktop\\A")
ann = read.table("circBank_circrna_annotation.txt", header  = T, sep = "\t")
ann$position = substr(ann$position, start = 4, stop = nchar(ann$position))
# install.packages("stringr")
library(stringr)
str1 = str_split(ann$position, ":", simplify = T) # simplify = T结果返回矩阵
str2 = str_split(str1[,2], "-", simplify = T)
# install.packages("magrittr")
library(magrittr)
str2[,1] = str2[,1] %>% as.numeric %>% +1
ann$position =  paste(str1[,1], str2[,1], sep = ":")
ann$position =  paste(ann$position, str2[,2], sep = "|")
ann = ann[,c(1,3)]
# write.table(ann, "circBankID.txt", sep = "\t")
data$circBankID = ann[match(data[,1], ann$position), 1]
data[,2] = gsub("n/a", "-", as.vector(data[,2]))
data[,4] = replace(data[,4], is.na(data[,4]), "-")
# ------------------------------------
#             三、合并数据
# ------------------------------------
dataname = grep("^data_", ls(), value = TRUE)
for(i in 1:length(dataname)) {
  temp = as.matrix(get(dataname[i])[,c(1,5)])
  data = cbind(data, temp[match(data[,1], temp[,1]), 2])
}
col_name = str_split(filename, "\\.", simplify = T)[,1]
colnum = ncol(data)
colnames(data)[5:colnum] = col_name
data = replace(data, is.na(data), 0)
keep = !duplicated(data[,1])
data = data[keep,]
write.table(data, "circRNA_count.txt", sep = "\t")
# ------------------------------------
#            四、RPM标化
# ------------------------------------
data[,5:colnum] = apply(data[,5:colnum], 2, as.numeric)
data[,5:colnum] = apply(data[,5:colnum], 2, function(x){exp(log(x) + log(1e6) - log(sum(x)))})
data = replace(data, is.na(data), 0)
# ------------------------------------
#            五、导出数据
# ------------------------------------
write.table(data, "circRNA_RPM.txt", sep = "\t")

