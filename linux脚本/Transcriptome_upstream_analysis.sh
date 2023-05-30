#!/bin/bash
# ------------------------------------
#              Trimmomatic
# ------------------------------------
### 1、单端测序
# case $1 in
# "1")
# echo -e " ------开始进行Trimmomatic !!!------ \n\n\n"
# mkdir -p ./result/Trimmomatic/$2
# Num=`ls ./rawdata/FTP/$2/*fastq.gz | wc -l` # 统计要处理的样本数，用以后面显示目前还剩几个要处理
# ls ./rawdata/FTP/$2/*fastq.gz | while read id ; do 
# 	trimmomatic SE -threads 15 -phred33 \
# 		$id \
# 		./result/Trimmomatic/$2/$(echo $id | awk -F '/' '{print $5}' | awk -F '.' '{print $1}').fq.gz \
# 		ILLUMINACLIP:/home/wucc/anaconda3/envs/RNA-seq/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa:2:30:10:8:FALSE \
# 		LEADING:3 TRAILING:3 \
# 		SLIDINGWINDOW:4:15 \
# 		MINLEN:$3
		# TruSeq3-SE.fa:污染序列文件，TruSeq2(GAII),TruSeq3(HiSeq MiSeq)
		# 2:seed搜索时允许的错配碱基个数
		# 30:PE palindrome模式下，R1R2之间进行接头切除所需比对分值
		# 10:切除接头序列最低比对分值，通常7-15
		# 8:PE palindrome模式下可切除的接头序列最短长度，历史原因默认8，实际上可以切至1bp的接头污染,2为官方设置
		# FALSE:整条去除与R1完全反向互补的R2,TRUE相反,True比较有用,会保留冗余信息，但是让流程更易于管理
		# LEADING:3 TRAILING:3:reads起始和末端端开始切除质量值低于3的碱基直到有个碱基其质量值达到阈值
		# SLIDINGWINDOW:4:15:滑窗剪切，从read5'端扫描切点，设置4bp窗口，碱基平均质量低于15就被切除
		# MINLEN:36:去掉处理后短于36的序列,一般读段都在100bp左右,默认36就好。如果读段是50bp甚至更短需要修改，改的越低，结果里就有越多的错误读段
# 	Num=$[$Num-1]
# 	echo -e "\n ------还剩 $Num 个样本 !!!------ \n"
# done
# echo -e "\n\n ------Trimmomatic结束 !!!------ \n\n\n"
# ;;
### 2、双端测序
# "2")
# echo -e " ------开始进行Trimmomatic !!!------ \n\n\n"
# mkdir -p ./result/Trimmomatic/$2
# Num=$[`ls ./rawdata/FTP/$2/*fastq.gz | wc -l` / 2] # 两个一起处理，所以文件数除以2，得到要处理的次数，用以后面显示目前还剩几个要处理
# ls ./rawdata/FTP/$2/*1.fastq.gz | while read id ; do
# 	Name=$(echo $id | awk -F '/' '{print $5}' | awk -F '_' '{print $1}') # 取出样本名除后缀外的内容
# 	trimmomatic PE -threads 15 -phred33 \
# 		./rawdata/FTP/$2/"$Name"_1.fastq.gz ./rawdata/FTP/$2/"$Name"_2.fastq.gz \
#         	-baseout ./result/Trimmomatic/$2/$Name.fq.gz \
#         	ILLUMINACLIP:/home/wucc/anaconda3/envs/RNA-seq/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10:2:TRUE \
#         	LEADING:3 TRAILING:3 \
#         	SLIDINGWINDOW:4:15 \
#         	MINLEN:$3
 		# -baseout:指定输出文件的basename，软件会自动为四个输出文件命名
 		# 输出文件：两个为paired（包含一对均经过筛选后留下的读段）和两个unpaired（包含仅有一个读段成功经过筛选）, 一般情况下若paired百分比占90%以上可只对其进行比对分析
# 	Num=$[$Num-1]
#         echo -e "\n ------还剩 $Num 个样本 !!!------ \n"
# done
# rm ./result/Trimmomatic/$2/*U.fq.gz
# echo -e "\n\n ------Trimmomatic结束 !!!------ \n\n\n"
# ;;
# esac
# ------------------------------------
#              一、fastp
# ------------------------------------
### 1、单端测序
case $1 in
"1")
echo -e "\n\n\n\e[32m\e[1m------fastp------\e[0m\n\n\n"
mkdir -p ./result/fastp/$2
Num=`ls ./rawdata/FTP/$2/*fastq.gz | wc -l` # 统计要处理的样本数，用以后面显示目前还剩几个要处理
for i in `ls ./rawdata/FTP/$2/*fastq.gz` ; do
	i=${i/.fastq.gz/}
	fastp -w 15 \
        	-i ${i}.fastq.gz -o ./result/fastp/$2/${i##*/}.fq.gz \
		-f $3 \
        	-j ./result/fastp/$2/${i##*/}.json -h ./result/fastp/$2/${i##*/}.html
	Num=$[$Num-1]
	echo -e "\n\e[32m\e[1m进度：还剩 $Num 个样本\e[0m\n"
done
echo -e "\n\n\e[32m\e[1mEnd：fastp结束\e[0m\n\n\n"
;;
### 2、双端测序
"2")
echo -e "\n\n\n\e[32m\e[1m------fastp------\e[0m\n\n\n"
mkdir -p ./result/fastp/$2
Num=`ls ./rawdata/FTP/$2/*1.fastq.gz | wc -l`
for i in `ls ./rawdata/FTP/$2/*1.fastq.gz` ; do
	i=${i/_1.fastq.gz/}
        fastp -w 15 \
        	-i ${i}_1.fastq.gz -o ./result/fastp/$2/${i##*/}_1.fq.gz \
		-I ${i}_2.fastq.gz -O ./result/fastp/$2/${i##*/}_2.fq.gz \
		-f $3 -F $3 \
        	-j ./result/fastp/$2/${i##*/}.json -h ./result/fastp/$2/${i##*/}.html
	Num=$[$Num-1]
        echo -e "\n\e[32m\e[1m进度：还剩 $Num 个样本\e[0m\n"
done
echo -e "\n\n\e[32m\e[1mEnd：fastp结束\e[0m\n\n\n"
;;
esac
# ------------------------------------
#           二、fastp结果统计
# ------------------------------------
echo -e "\e[32m\e[1m------统计fastp结果------\e[0m\n\n\n"
Rscript fastp_json.R -p "~/Transcriptome_upstream_analysis/result/fastp/$2/"
for i in `cat ./result/fastp/$2/del_sample.txt` ; do
	rm ./result/fastp/$2/$i*
done
echo -e "\n\n\n\e[32m\e[1mEnd：fastp结果统计完成\e[0m\n\n\n"
# ------------------------------------
#          三、STAR--建立索引
# ------------------------------------
if [ ! -d ./result/STAR/index/$4 ] ; then
	echo -e "\e[32m\e[1m------建立 readlength=$4 的索引------\e[0m\n\n\n"
	mkdir -p ./result/STAR/index/$4
	len=$[$4-1]
	STAR \
        	--runThreadN 15 \
        	--runMode genomeGenerate \
        	--genomeDir ./result/STAR/index/$4 \
        	--genomeFastaFiles ./rawdata/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
        	--sjdbGTFfile ./rawdata/GRCh37/Homo_sapiens.GRCh37.75.gtf \
        	--sjdbOverhang $len  
		# --runMode:工作模式，构建基因组索引
		# --genomeDir:构建好的索引所在目录，需提前建好，绝对路径比较好，接下来也是
		# --genomeFastaFiles、--sjdbGTFfile：参考基因组和注释信息位置，要同一个版本
		# --sjdbOverhang；默认值是100，理想值read长度减1，junction周围的序列长度
	echo -e "\n\n\n\e[32m\e[1mEnd：readlength=$4 的索引建立完成\e[0m\n\n\n"
else
	echo -e "\e[32m\e[1m------readlength=$4 的索引已经建立------\e[0m\n\n\n"
fi
# ------------------------------------
#            四、STAR--比对
# ------------------------------------
### 1、单端
case $1 in
"1")
echo -e "\e[32m\e[1m------STAR------\e[0m\n\n\n"
mkdir -p ./result/STAR/align/$2
len=$[$4-1]
Num=`ls ./result/fastp/$2/*.fq.gz | wc -l`
for i in `ls ./result/fastp/$2/*.fq.gz` ; do
	i=${i/.fq.gz/}
    	STAR \
    		--runThreadN 15 \
                --runMode alignReads \
                --genomeDir ./result/STAR/index/$4 \
                --sjdbGTFfile ./rawdata/GRCh37/Homo_sapiens.GRCh37.75.gtf \
		--readFilesIn ${i}.fq.gz \
		--outFileNamePrefix ./result/STAR/align/$2/${i##*/}_ \
		--readFilesCommand zcat \
    		--outSAMtype BAM Unsorted \
		--outReadsUnmapped Fastx \
		--quantMode GeneCounts \
                --sjdbOverhang $len 
		# --readFilesCommand zcat：对所有gz文件解压
                # --outSAMtype BAM Unsorted: 将SAM格式转换为BAM格式，不排序，实际上就是按name的sort，下游可以直接接HTSeq；默认SAM；第一个取值SAM或BAM，第二个值Unsorted或SortedByCoordinate
		# --outReadsUnmapped Fastx：输出未比对上的reads
                # --quantMode GeneCounts：计数每个基因的reads数，需要--sjdbGTFfile指定GTF文件；中间加TranscriptomeSAM将基因组比对结果转换为转录组结果	
	Num=$[$Num-1]
        echo -e "\n\e[32m\e[1m进度：还剩 $Num 次比对\e[0m\n"
done
rm ./result/STAR/align/$2/{*.out,*SJ*}
rm -rf ./result/STAR/align/$2/*genome
echo -e "\n\n\e[32m\e[1mEnd：STAR比对结束\e[0m\n\n\n"
;;  
### 2、双端
"2")
echo -e "\e[32m\e[1m------STAR------\e[0m\n\n\n"
mkdir -p ./result/STAR/align/$2
len=$[$4-1] 
Num=`ls ./result/fastp/$2/*1.fq.gz | wc -l`
for i in `ls ./result/fastp/$2/*1.fq.gz` ; do
	i=${i/_1.fq.gz/}
    	STAR \
        	--runThreadN 15 \
        	--runMode alignReads \
        	--genomeDir ./result/STAR/index/$4 \
        	--sjdbGTFfile ./rawdata/GRCh37/Homo_sapiens.GRCh37.75.gtf \
		--readFilesIn ${i}_1.fq.gz ${i}_2.fq.gz \
		--outFileNamePrefix ./result/STAR/align/$2/${i##*/}_ \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
		--outReadsUnmapped Fastx \
		--quantMode GeneCounts \
		--sjdbOverhang $len
        Num=$[$Num-1]
        echo -e "\n\e[32m\e[1m进度：还剩 $Num 次比对\e[0m\n"
done
rm ./result/STAR/align/$2/{*.out,*SJ*}
rm -rf ./result/STAR/align/$2/*genome
echo -e "\n\n\e[32m\e[1mEnd：STAR比对结束\e[0m\n\n\n"
;;
esac
# ------------------------------------
#              htseq-count
# ------------------------------------
### 1、bam文件按name排序
# echo -e " ------开始将bam文件按name排序 !!!------ \n\n\n"
# for f in `ls ./result/STAR/align/$2/*_Aligned.out.bam`; do mv $f ${f%_*}.Bam ; done # 方便之后并行时命名输出文件
# ls ./result/STAR/align/$2/*.Bam | parallel --eta -j 60% 'samtools sort -n {} -o {.}.sorted' # --eta显示任务完成的预计剩余时间;-j 60%使用60%的核心数
# rm ./result/STAR/align/$2/*.Bam
# echo -e "\n\n\n ------bam文件按name排序完成 !!!------ \n\n\n"
### 2、计数
# echo -e " ------开始进行htseq-count计数 !!!------ \n\n\n"
# mkdir -p ./result/htseq_count/$2
# ls ./result/STAR/align/$2/*.sorted | parallel --eta -j 60% 'htseq-count -f bam -r name -s no -q {} ./rawdata/GRCh38/Homo_sapiens.GRCh38.108.gtf > ./result/htseq_count/$2/{/.}.count' # -s是否这个数据是来自链特异性建库,默认yes
# echo -e "\n\n\n ------htseq-count计数完成 !!!------ \n\n\n"
# ------------------------------------
#            五、合并矩阵
# ------------------------------------
echo -e "\e[32m\e[1m------合并矩阵------\e[0m\n\n\n"
mkdir -p ./result/matrix/$2
echo "gene_id" |  cat - `ls ./result/STAR/align/$2/*ReadsPerGene* | head -n 1` | sed "2,5d" | colrm 16 > ./result/matrix/$2/name.count # colrm 删除列，这边指删除第16列以后的所有；此行为了取出行名
ls ./result/STAR/align/$2/*ReadsPerGene* | while read id ; do
	pre=`basename $id _ReadsPerGene.out.tab` # 取出文件名，不包括后缀
	cat $id | sed "1,4d" | awk -F ' ' '{print $2}' > ./result/matrix/$2/mid.count
	echo $pre | cat -  ./result/matrix/$2/mid.count > ./result/matrix/$2/${pre}_count
done
paste ./result/matrix/$2/name.count ./result/matrix/$2/*_count > ./result/matrix/$2/rawdata
rm ./result/matrix/$2/*count
echo -e "\e[32m\e[1mEnd：合并完成\e[0m\n\n\n"
# ------------------------------------
#            六、rawdata处理
# ------------------------------------
echo -e "\e[32m\e[1m------处理rawdata------\e[0m\n\n\n"
Rscript Standardize_classify.R -p "~/Transcriptome_upstream_analysis/result/matrix/$2/" -a "./rawdata/GRCh37/"
echo -e "\n\n\n\e[32m\e[1mEnd：处理完成\e[0m\n\n\n"
# ------------------------------------
#          七、构建BWA索引
# ------------------------------------
# mkdir -p ./result/BWA/index
# cd ./result/BWA/index
# bwa index \
#       -p hg19_index \
#       -a bwtsw \
#       ~/Transcriptome_upstream_analysis/rawdata/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
#       # -p输出数据库的前缀; -a构建index的算法
# cd ~/Transcriptome_upstream_analysis/
# ------------------------------------
#             八、BWA比对
# ------------------------------------
### 1、单端测序
case $1 in
"1")
echo -e "\e[32m\e[1m------BWA比对------\e[0m\n\n\n"
Num=`ls ./result/STAR/align/$2/*Unmapped* | wc -l`
mkdir -p ./result/BWA/align/$2
for i in `ls ./result/STAR/align/$2/*Unmapped*` ; do
        Num=$[$Num-1]
        pre=`basename $i _Unmapped.out.mate1` # 取出文件名，不包括后缀
        bwa mem \
                -T 19 \
                -t 15 \
                ./result/BWA/index/hg19_index \
                $i \
                > ./result/BWA/align/$2/$pre.sam
        echo -e "\n\e[32m\e[1m进度：还剩 $Num 个样本\e[0m\n"
done
echo -e "\n\n\e[32m\e[1mEnd：比对结束\e[0m\n\n\n"
;;
### 2、双端测序
"2")
echo -e "\e[32m\e[1m------BWA比对------\e[0m\n\n\n"
Num=`ls ./result/STAR/align/$2/*Unmapped.out.mate1 | wc -l`
mkdir -p ./result/BWA/align/$2
for i in `ls ./result/STAR/align/$2/*Unmapped.out.mate1` ; do
        Num=$[$Num-1]
        pre=`basename $i _Unmapped.out.mate1`
        bwa mem \
                -T 19 \
                -t 15 \
                ./result/BWA/index/hg19_index \
                $i ${i%.*}.mate2 \
                > ./result/BWA/align/$2/$pre.sam
        echo -e "\n\e[32m\e[1m进度：还剩 $Num 个样本\e[0m\n"
done
echo -e "\n\n\e[32m\e[1mEnd：比对结束\e[0m\n\n\n"
;;
esac
# ------------------------------------
#              九、CIRI2
# ------------------------------------
echo -e "\e[32m\e[1m------CIRI2------\e[0m\n\n\n"
Num=`ls ./result/BWA/align/$2/*sam | wc -l`
mkdir -p ./result/CIRI2/$2
for i in `ls ./result/BWA/align/$2/*sam` ; do
        Num=$[$Num-1]
        pre=`basename $i .sam`
        perl ./tool/CIRI/CIRI_v2.0.6/CIRI2.pl \
                -T 10 \
                -F ./rawdata/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
                -A ./rawdata/GRCh37/Homo_sapiens.GRCh37.75.gtf \
                -I $i \
                -O ./result/CIRI2/$2/$pre.txt \
                -Q
        echo -e "\n\e[32m\e[1m进度：还剩 $Num 个样本\e[0m\n"
done
rm ./result/CIRI2/$2/*.log
echo -e "\n\n\e[32m\e[1mEnd：CIRI2结束\e[0m\n\n\n"
# ------------------------------------
#          十、circRNA结果处理
# ------------------------------------
echo -e "\e[32m\e[1m------处理circRNA结果------\e[0m\n\n\n"
mkdir -p ./result/circRNA/$2
Rscript circRNA.R -p "./result/CIRI2/$2/" -a "~/Transcriptome_upstream_analysis/rawdata/circBankID/" -o "~/Transcriptome_upstream_analysis/result/circRNA/$2"
echo -e "\n\n\n\e[32m\e[1mEnd：处理完成\e[0m\n\n\n"
