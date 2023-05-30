# ------------------------------------
#           一、导入gtf文件
# ------------------------------------
setwd("C:\\Users\\13395\\Desktop")
# install.packages("rtracklayer")
library(rtracklayer)
gtf = rtracklayer::import("Homo_sapiens.GRCh37.75.gtf")
gtf = as.data.frame(gtf)
# write.table(gtf, "gtf.txt", sep = "\t")
# gtf = read.table("gtf.txt", header = T, sep = "\t") 
# ------------------------------------
#            二、导入count
# ------------------------------------
count = read.table("rawdata", header = T, sep = "\t") 
rownames(count) = count[,1]
count = count[,-1]
# ------------------------------------
#         三、gtf文件预处理
# ------------------------------------
exon = gtf[gtf$type == "exon", c("start", "end", "gene_id")]
gl = lapply(split(exon, exon$gene_id), function(x) { # split按gene_id对exon进行分割成list；lapply按list进行function操作
  tmp=apply(x, 1, function(y){y[1]:y[2]}) # apply对每个list按行进行function
  length(unique(unlist(tmp)))
  }) # 计算各基因的外显子去冗余加和
gl1 = data.frame(gene_id = names(gl), length = as.numeric(gl))
inter = intersect(rownames(count), gl1$gene_id)
gl1 = gl1[match(inter, gl1$gene_id),]
# ------------------------------------
#             四、TPM标化
# ------------------------------------
tmp = apply(count, 2, function(x){x / gl1$length})
TPM = apply(tmp, 2, function(x){exp(log(x) + log(1e6) - log(sum(x)))})
TPM = cbind(rownames(TPM), TPM)
colnames(TPM)[1]="gene_id"
# ------------------------------------
#            五、分类去重
# ------------------------------------
### 1、mRNA
tmp = gtf[which(gtf$type == "gene"), c("gene_id", "gene_name", "gene_biotype")]
mRNA_id = as.matrix(tmp[which(tmp$gene_biotype == "protein_coding"), c("gene_id", "gene_name")])
mRNA = merge(mRNA_id, TPM, by = "gene_id")
mRNA_count = count[match(mRNA[,1], rownames(count)),]
mRNA = mRNA[,-1]
mRNA[,-1] = apply(mRNA[,-1], 2, as.numeric) # 按列转变为数值型
index = order(rowMeans(mRNA[,-1]), decreasing = T) # 计算行平均值，按降序排列
mRNA = mRNA[index,] # 调整表达谱的基因顺序
keep = !duplicated(mRNA[,1]) # 对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
mRNA = mRNA[keep,] # 得到处理之后的表达谱矩阵
rownames(mRNA) = mRNA[,1]
mRNA_TPM = mRNA[,-1]
if (length(which(apply(mRNA_count, 1, sum) == 0)) != 0) {
  mRNA_count = mRNA_count[-which(apply(mRNA_count, 1, sum) == 0),]
}
if (length(which(apply(mRNA_TPM, 1, sum) == 0)) != 0) {
  mRNA_TPM = mRNA_TPM[-which(apply(mRNA_TPM, 1, sum) == 0),]
}
### 2、lncRNA
lncRNA_id = as.matrix(tmp[which(tmp$gene_biotype == "3prime_overlapping_ncrna" | tmp$gene_biotype == "antisense" | tmp$gene_biotype == "sense_intronic" | tmp$gene_biotype == "sense_overlapping" | tmp$gene_biotype == "lincRNA"), c("gene_id", "gene_name")])
lncRNA = merge(lncRNA_id, TPM, by="gene_id")
lncRNA_count = count[match(lncRNA[,1], rownames(count)),]
lncRNA = lncRNA[,-1]
lncRNA[,-1] = apply(lncRNA[,-1], 2, as.numeric) 
index = order(rowMeans(lncRNA[,-1]), decreasing = T) 
lncRNA = lncRNA[index,]
keep = !duplicated(lncRNA[,1]) 
lncRNA = lncRNA[keep,] 
rownames(lncRNA) = lncRNA[,1]
lncRNA_TPM = lncRNA[,-1
if (length(which(apply(lncRNA_count, 1, sum) == 0)) != 0) {
  lncRNA_count  = lncRNA_count[-which(apply(lncRNA_count , 1, sum) == 0),]
}
if (length(which(apply(lncRNA_TPM, 1, sum) == 0)) != 0) {
  lncRNA_TPM  = lncRNA_TPM[-which(apply(lncRNA_TPM , 1, sum) == 0),]
}
# ------------------------------------
#            六、导出数据
# ------------------------------------
write.table(mRNA_count, "mRNA_count.txt", sep = "\t")
write.table(lncRNA_count, "lncRNA_count.txt", sep = "\t")
write.table(mRNA_TPM, "mRNA_TPM.txt", sep = "\t")
write.table(lncRNA_TPM, "lncRNA_TPM.txt", sep = "\t")

