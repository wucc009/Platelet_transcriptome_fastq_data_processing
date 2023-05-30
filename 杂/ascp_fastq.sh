#!/bin/bash
#运行开始及时间
echo -e "\n\n\n ascp download fastq begin !!! \n\n\n"
date
#创建切换输出储存目录
mkdir -p ~/Test/rawdata/fq/
cd ~/Test/rawdata/fq/
pwd
#创建变量
openssh=~/.aspera/connect/etc/asperaweb_id_dsa.openssh #密钥路径
#获取数据
cat ~/Test/SRR_Acc_List.txt | while read id
do
num=`echo $id | wc -m ` #wc -m 会把每行结尾$也算为一个字符
echo "SRRnum+1 is $num"
#样本编号为SRR+8位数
if [ $num -eq 12 ]
then
        echo "SRR + 8"
        x=$(echo $id | cut -b 1-6)
        y=$(echo $id | cut -b 10-11)
        echo "Downloading $id "
        ( ascp  -QT -l 500m -P33001  -k 1 -i $openssh \
                era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/0$y/$id/   ./ & )
#如果样本编号为SRR+7位数 #
elif [ $num -eq 11 ]
then
        echo  "SRR + 7"
        x=$(echo $id | cut -b 1-6)
        y=$(echo $id | cut -b 10-10)
        echo "Downloading $id "
        ( ascp  -QT -l 500m -P33001  -k 1 -i $openssh \
                        era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/00$y/$id/   ./ & )
#如果样本编号为SRR+6位数 #
elif [ $num -eq 10 ]
then
        echo  "SRR + 6"
        x=$(echo $id |cut -b 1-6)
        echo "Downloading $id "
        ( ascp  -QT -l 500m -P33001 -k 1 -i  $openssh \
                        era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/$id/   ./ & )
fi
done
