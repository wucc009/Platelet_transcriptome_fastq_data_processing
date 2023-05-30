#!/home/wucc/anaconda3/envs/RNA-seq/bin/R
# ------------------------------------
#               设置参数
# ------------------------------------  
cat("\n\033[32m\033[1mStep 1: 设置参数\033[0m")
# install.packages('optparse')
require('optparse')
option_list = list( # 构建参数列表	
	make_option(c("-p", "--path"), type = "character", default = FALSE, help = "输入引号括起来的rawdata路径"))
	# type为输入参数的类型，有“logical”, “integer”, “double”, “complex”, “character”这几种类型 
opt = parse_args(OptionParser(option_list = option_list)) # 解析参数 
# ------------------------------------
#               统计结果
# ------------------------------------  
cat("\n\033[32m\033[1mStep 2: 统计结果\033[0m")
setwd(opt$path)
require(jsonlite)
file_name = c(list.files(pattern = "json$")) 
data = NULL
for(i in 1:length(file_name)) { 
  jsonData = fromJSON(file_name[i])
  temp = cbind(unlist(strsplit(file_name[i], "\\."))[1], jsonData[[1]][[3]][[1]], jsonData[[1]][[4]][[1]], jsonData[[1]][[4]][[6]])
  data = rbind(data, temp)
}
colnames(data) = c("Sample_ID", "raw_reads", "clean_reads", "%_>Q30")
del_sample = c(data[which(data[,4] < 0.85), 1])
write.table(data, "json_result.txt", sep = "\t")
write.table(del_sample, "del_sample.txt", sep = "\t", quote = F, row.names = F, col.names = F)
