#!/bin/bash
# ------------------------------------
#                二合一
# ------------------------------------
case $1 in
"2")
echo -e "\n\n\n ------开始进行数据二合一 !!!------ \n\n\n"
Num=`wc -l < ./rawdata/FTP/merge_data/$2`
sed -e 's/.$//' ./rawdata/FTP/merge_data/$2 | cat | while read id ; do # sed -e可对一行进行多个操作；s代表替换substitute；.$代表行尾；s/.$//把行尾替换为空
        cat ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $2}'`.fastq.gz >> ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $1}'`.fastq.gz
        rm ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $2}'`.fastq.gz
	Num=$[$Num-1]
        echo -e " ------还剩 $Num 次合并 !!!------ \n"
done
echo -e "\n\n ------数据二合一结束 !!!------ \n\n\n"
;;
# ------------------------------------
#                四合一
# ------------------------------------
"4")
echo -e "\n\n\n ------开始进行数据四合一 !!!------ \n\n\n"
Num=`wc -l < ./rawdata/FTP/merge_data/$2`
sed -e 's/.$//' ./rawdata/FTP/merge_data/$2 | cat | while read id ; do
	cat ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $4}'`.fastq.gz >> ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $3}'`.fastq.gz
	cat ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $3}'`.fastq.gz >> ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $2}'`.fastq.gz
	cat ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $2}'`.fastq.gz >> ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $1}'`.fastq.gz
	rm ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $2}'`.fastq.gz ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $3}'`.fastq.gz ./rawdata/FTP/$3/`echo $id | awk -F " " '{print $4}'`.fastq.gz
	Num=$[$Num-1]
        echo -e " ------还剩 $Num 次合并 !!!------ \n"
done
echo -e "\n\n ------数据四合一结束 !!!------ \n\n\n"
;;
esac
