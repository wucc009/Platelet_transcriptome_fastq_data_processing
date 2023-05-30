#!/home/wucc/anaconda3/envs/RNA-seq/bin/R
# ------------------------------------
#             一、设置参数
# ------------------------------------
cat("\033[32m\033[1mStep 1: 设置参数\033[0m") # \033[32m绿色；\033[1m加粗；\033[0m复原 
# install.packages('optparse')
require('optparse')
option_list = list(  # 构建参数列表	
	make_option( c("-p", "--path"), type = "character", default = FALSE, help = "输入引号括起来的rawdata路径"),
	make_option( c("-a", "--annotation"), type = "character", default = FALSE, help = "输入引号括起来的注释文件路径"),
	make_option( c("-o", "--outpath"), type = "character", default = FALSE, help = "输入引号括起来的输出路径"))  
	# type为输入参数的类型，有“logical”, “integer”, “double”, “complex”, “character”这几种类型 
opt = parse_args(OptionParser(option_list = option_list)) # 解析参数 
# ------------------------------------
#       二、导入数据及预处理
# ------------------------------------
cat("\n\033[32m\033[1mStep 2: 导入数据及预处理\033[0m")
setwd(opt$path)
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
#          三、添加circBankID
# ------------------------------------
cat("\n\033[32m\033[1mStep 3: 添加circBankID\033[0m")
setwd(opt$annotation)
ann = read.table("circBankID.txt", header  = T, sep = "\t")
data$circBankID = ann[match(data[,1], ann$position), 1]
data[,2] = gsub("n/a", "-", as.vector(data[,2]))
data[,4] = replace(data[,4], is.na(data[,4]), "-")
# ------------------------------------
#             四、合并数据
# ------------------------------------
cat("\n\033[32m\033[1mStep 4: 合并数据\033[0m")
dataname = grep("^data_", ls(), value = TRUE)
for(i in 1:length(dataname)) {
  temp = as.matrix(get(dataname[i])[,c(1,5)])
  data = cbind(data, temp[match(data[,1], temp[,1]), 2])
}
require('stringr')
col_name = str_split(filename, "\\.", simplify = T)[,1]
colnum = ncol(data)
colnames(data)[5:colnum] = col_name
data = replace(data, is.na(data), 0)
keep = !duplicated(data[,1])
data = data[keep,]
setwd(opt$outpath)
write.table(data, "circRNA_count.txt", sep = "\t")
# ------------------------------------
#            五、RPM标化
# ------------------------------------
cat("\n\033[32m\033[1mStep 5: RPM标化\033[0m")
data[,5:colnum] = apply(data[,5:colnum], 2, as.numeric)
data[,5:colnum] = apply(data[,5:colnum], 2, function(x){exp(log(x) + log(1e6) - log(sum(x)))})
data = replace(data, is.na(data), 0)
# ------------------------------------
#            六、导出数据
# ------------------------------------
cat("\n\033[32m\033[1mStep 6: 导出数据\033[0m")
write.table(data, "circRNA_RPM.txt", sep = "\t")
