#/bin/bash
# ------------------------------------
#               单端测序
# ------------------------------------
### 1、判断单双端
case $1 in
"1")
### 2、产生临时文件
mkdir $3
less -S ./rawdata/FTP/$2 | sed '1d' > ./$3/sample.txt # -S一行展示;将输入文件去除第一行标签
awk '{print "ftp://"$3}' ./$3/sample.txt > ./$3/sample_ftp.txt # 取出FTP
awk -F '/' '{print $9}' ./$3/sample_ftp.txt > ./$3/sample_name.txt # 从FTP中取出文件名
awk '{print $2}' ./$3/sample.txt > ./$3/sample_md5.txt # 取出md5码
### 3、下载数据
echo -e "\n\n\n\e[32m\e[1m------下载数据------\e[0m\n\n\n" # -e将转义后的内容输出到屏幕上
mkdir ./rawdata/FTP/${2%.*} # ${2%.*}删除$2最后一个.及其右边所有内容
aria2c -c --max-download-limit=15M -d ./rawdata/FTP/${2%.*} -i ./$3/sample_ftp.txt # -c断点传续；--max-download-limit=15M限制上传带宽；-d设置目录；-i按文件作为输入
### 4、校验所得文件的md5码
echo -e "\n\n\n\e[32m\e[1m------校验md5值------\e[0m\n\n\n"
for line in $(cat ./$3/sample_name.txt) ; do echo ./rawdata/FTP/${2%.*}/$line >> ./$3/loc.txt ; done # 建立各样本储存路径
paste -d " " ./$3/sample_md5.txt ./$3/loc.txt > ./$3/md5.txt # -d " "以空格作为分隔符；建立MD5值及对应文件路径的文件
md5sum -c ./$3/md5.txt | grep -n "FAILED" > ./$3/diff.txt # -n输出行号；检查文件的MD5值,并找出差异
### 5、若存在差异，则进行以下循环直到全部下完整为止，无差异则跳过循环下载结束
NUM=1 # 区别第几次循环
while [ `head -1 ./$3/diff.txt` ] # 只看第一行是因为无差异为空，不进行以下循环;有差异的话会得到格式奇怪的结果,while [ `cat ./$3/diff.txt` ]会报错，取第一行判断不为空即进行以下循环
do
	echo -e "\n\n\n\e[32m\e[1m------重新下载不完整数据 第 $NUM 次------\e[0m\n\n\n"
	cat ./$3/diff.txt | awk -F ':' '{print $1}' > ./$3/order.txt # 以:分割取第一部分得到差异行号
	for i in $(cat ./$3/order.txt); do awk -v a=$i 'NR==a' ./$3/sample_name.txt >> ./$3/sample_name_1.txt; done # 取出sample_name.txt中的这些行
	rm `awk '{print "./rawdata/FTP/${2%.*}/"$1}' ./$3/sample_name_1.txt` # 删除下载不完整的数据
	for i in $(cat ./$3/order.txt); do awk -v a=$i 'NR==a' ./$3/sample_ftp.txt >> ./$3/sample_ftp_1.txt; done # -v创建变量;取出sample_ftp中的这些行	
	aria2c -c --max-download-limit=15M -d ./rawdata/FTP/${2%.*} -i ./$3/sample_ftp_1.txt # 重新下载
	echo -e "\n\n\n\e[32m\e[1m------重新校验md5值 第 $NUM 次------\e[0m\n\n\n"
        for i in $(cat ./$3/order.txt); do awk -v a=$i 'NR==a' ./$3/md5.txt >> ./$3/md5_1.txt; done # 取出md5.txt中的这些行
	md5sum -c ./$3/md5_1.txt | grep -n "FAILED" > ./$3/diff.txt	
	cat ./$3/sample_name_1.txt > ./$3/sample_name.txt # 替换数据
	cat ./$3/sample_ftp_1.txt > ./$3/sample_ftp.txt
	cat ./$3/md5_1.txt > ./$3/md5.txt
	rm ./$3/sample_name_1.txt ./$3/sample_ftp_1.txt ./$3/md5_1.txt # 因为这些文件是追加的，所以下一循环前得删除先
	NUM=$[$NUM+1]
done
echo -e "\e[32m\e[1m------无不完整数据，下载完成------\e[0m\n\n\n"
rm -rf ./$3
;;
# ------------------------------------
#               双端测序
# ------------------------------------
"2")
mkdir $3
less -S ./rawdata/FTP/$2 | sed '1d' > ./$3/sample.txt
awk '{print $3}' ./$3/sample.txt | awk -F ';' '{print "ftp://"$1"\n""ftp://"$2}' > ./$3/sample_ftp.txt
awk -F '/' '{print $9}' ./$3/sample_ftp.txt > ./$3/sample_name.txt
awk '{print $2}' ./$3/sample.txt | awk -F ';' '{print $1"\n"$2}' > ./$3/sample_md5.txt
echo -e "\n\n\n\e[32m\e[1m------下载数据------\e[0m\n\n\n"
mkdir ./rawdata/FTP/${2%.*}
aria2c -c --max-download-limit=15M -d ./rawdata/FTP/${2%.*} -i ./$3/sample_ftp.txt
echo -e "\n\n\n\e[32m\e[1m------校验md5值------\e[0m\n\n\n"
for line in $(cat ./$3/sample_name.txt) ; do echo ./rawdata/FTP/${2%.*}/$line >> ./$3/loc.txt ; done
paste -d " " ./$3/sample_md5.txt ./$3/loc.txt > ./$3/md5.txt
md5sum -c ./$3/md5.txt | grep -n "FAILED" > ./$3/diff.txt
NUM=1
while [ `head -1 ./$3/diff.txt` ]
do
	echo -e "\n\n\n\e[32m\e[1m------重新下载不完整数据 第 $NUM 次------\e[0m\n\n\n"
	cat ./$3/diff.txt | awk -F ':' '{print $1}' > ./$3/order.txt
	for i in $(cat ./$3/order.txt); do awk -v a=$i 'NR==a' ./$3/sample_name.txt >> ./$3/sample_name_1.txt; done
	rm `awk '{print "./Transcriptome_upstream_analysis/"$1}' ./$3/sample_name_1.txt`	
	for i in $(cat ./$3/order.txt); do awk -v a=$i 'NR==a' ./$3/sample_ftp.txt >> ./$3/sample_ftp_1.txt; done
	aria2c -c --max-download-limit=15M -d ./rawdata/FTP/${2%.*} -i ./$3/sample_ftp_1.txt
	echo -e "\n\n\n\e[32m\e[1m------重新校验md5值 第 $NUM 次------\e[0m\n\n\n"
	for i in $(cat ./$3/order.txt); do awk -v a=$i 'NR==a' ./$3/sample_md5.txt >> ./$3/sample_md5_1.txt; done
	md5sum -c ./$3/md5_1.txt | grep -n "FAILED" > ./$3/diff.txt
	cat ./$3/sample_name_1.txt > ./$3/sample_name.txt  	
	cat ./$3/sample_ftp_1.txt > ./$3/sample_ftp.txt
	cat ./$3/md5_1.txt > ./$3/md5.txt
	rm ./$3/sample_md5_1.txt ./$3/sample_ftp_1.txt ./$3/md5_1.txt
	NUM=$[$NUM+1]
done
echo -e "\e[32m\e[1m------无不完整数据，下载完成------\e[0m\n\n\n"
rm -rf ./$3
;;
esac
