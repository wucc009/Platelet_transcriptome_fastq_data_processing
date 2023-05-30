setwd("C:\\Users\\13395\\Desktop")
# install.packages("jsonlite")
library(jsonlite)
file_name = c(list.files(pattern = "json$"))
data = NULL
for(i in 1:length(file_name)) {
  jsonData = fromJSON(file_name[i])
  temp = cbind(unlist(strsplit(file_name[i], "\\."))[1], jsonData[[1]][[3]][[1]], jsonData[[1]][[4]][[1]], jsonData[[1]][[4]][[6]])
  data = rbind(data, temp)
}
colnames(data) = c("Sample_ID", "raw_reads", "clean_reads", "%_>Q30")
del_sample = c(data[which(data[,4] < 0.85), 1])
write.table(data, "json_result.txt", sep="\t")
write.table(del_sample, "del_sample.txt", sep="\t", quote=F, row.names=F, col.names=F)
