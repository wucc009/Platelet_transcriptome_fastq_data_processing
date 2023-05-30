# Platelet_transcriptome_fastq_data_processing
这是一个用于转录组fastq数据的处理管道
## 特点
- 可进行单双端fastq数据的从头分析
- 最终获得mRNA、lncRNA、circRNA表达数据
- 只需少量操作即可完成
## 使用指南
1. 事前准备：
- 创建一个虚拟环境，安装所需软件及包（此处不一一列举，可在代码中获取）
- 从ensemble数据库下载所需.fa和.gtf文件
- 从circbank下载所需文件
2. 先从GEO、SRA数据库收集数据，最后从ENA获得.txt文件，包含三列：run_accession、fastq_md5、fastq_ftp：
  - [GEO](https://www.ncbi.nlm.nih.gov/geo/?tdsourcetag=s_pcqq_aiomsg)
  - [SRA](https://www.ncbi.nlm.nih.gov/sra)
  - [ENA](https://www.ebi.ac.uk/ena/browser/home)
---
3. 通过调用Aria2c_fastq.sh脚本进行下载：

`nohup bash Aria2c_fastq.sh 1 01-PRJNA722042-3-single.txt temp_file1 >Aria2c_log1 2>&1 &`
- 1：单端（2即双端）
- 01-PRJNA722042-3-single.txt：第一步下载的文件
- temp_file1：程序运行中的临时文件夹
---
4. 使用fastq、fastqc和fastp进行质控：
- fastq代码不小心删了（自己写一个很简单）
- 之后调用fastQC.sh

`nohup bash fastQC.sh 01-PRJNA722042-3-single.txt log1 2>&1 &`
- 根据fastQC结果最后使用Transcriptome_upstream_analysis.sh中的fastp进行质控
---
5. 调用Transcriptome_upstream_analysis.sh脚本：
- 进行STAR比对，获得rawdata
- 对rawdata合并，并进行TPM标化，然后输出mRNA、lncRNA表达数据
- 进行BWA比对
- 进行CIRI2获取circRNA表达数据

`nohup bash Transcriptome_upstream_analysis.sh 1 01-PRJNA722042-3-single 3 47 >log1 2>&1 &`
- 3:fastp对f端修剪长度
- 47:star中参数--sjdbOverhang的值，默认值是100，理想值为read长度减1，这边采用read长度减1
