export PERL5LIB=/Share/home/zhangqf7/gongjing/software/czplib

sample=$1
# cd /Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/${sample}
cd /150T/zhangqf7/jinsong_sequencing_rawdata/modification/20190623-iclip/Rawdata/${sample}

fq=${sample}_R1.fq.gz

cutadapt -a AGATCGGAAGAGCACACGTCTGAAC -o cutadapt-all.fastq.gz $fq
gunzip -c cutadapt-all.fastq.gz > cutadapt-all.fastq

perl /Share/home/zhangqf7/gongjing/software/ctk-master/fastq2collapse.pl cutadapt-all.fastq deduplicate_total.fastq 
rm cutadapt-all.fastq


#########################################################
##去tail A
##Strip barcode
##去barcode

cutadapt -a "A{10}" -e 0.1 -o rmPOLYA_total.fastq deduplicate_total.fastq
awk '{if(NR%4==1){split($1,a,"#");print a[1]"#"a[2]}else{print $0}}' rmPOLYA_total.fastq > temp1.fastq
gzip deduplicate_total.fastq

perl /Share/home/zhangqf7/gongjing/software/ctk-master/stripBarcode.pl -format fastq -len 3 temp1.fastq debarcode_total.fastq
rm temp1.fastq
gzip rmPOLYA_total.fastq


#########################################################
##Read quality filtering
##测序reads质量控制
##18nt保证mapping准确率

java -Xmx4g -jar /Share/home/zhangqf7/gongjing/zebrafish/script/icSHAPE-master/bin/trimmomatic-0.30.jar SE -phred33 debarcode_total.fastq qf_trim_total.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:18
gzip debarcode_total.fastq

###
mkdir fastqc
/Share/home/zhangqf/shaodi/app/fastqc/fastqc qf_trim_total.fastq -t 3 -o ./fastqc
/Share/home/zhangqf/shaodi/app/fastqc/fastqc rmPOLYA_total.fastq.gz -t 3 -o ./fastqc
/Share/home/zhangqf/shaodi/app/fastqc/fastqc debarcode_total.fastq.gz -t 3 -o ./fastqc

######mapping
genome_index=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/genome/bwa/danRer10.fa
gtf=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/danRerV10.refGene.gtf

mkdir bwa
cd bwa
/Share/home/zhangqf7/gongjing/software/bwa-0.7.17/bwa aln -t 20 -n 0.06 -q 20 ${genome_index} ../qf_trim_total.fastq > iCLIP.sai
/Share/home/zhangqf7/gongjing/software/bwa-0.7.17/bwa samse ${genome_index} iCLIP.sai ../qf_trim_total.fastq > iCLIP.sam
    
samtools view iCLIP.sam -q 20 -h -S > uniqmap.sam
echo "reads with no -q:" > readsCount
cat iCLIP.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
samtools view -S iCLIP.sam -b -o iCLIP.bam

echo "reads with uniqmap:" >> readsCount
cat uniqmap.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount

samtools view -S uniqmap.sam -b -o uniqmap.bam
bedtools bamtobed -i uniqmap.bam -split > uniqmap.bed

echo "" >> readsCount
echo "" >> readsCount
echo "" >> readsCount

echo "" >> readsCount
echo "" >> readsCount
echo "" >> readsCount

cd /Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/${sample}/bwa
CTK_N=CTK_Procedure1
mkdir -p ${CTK_N}
cd ${CTK_N}

#####1. Parsing SAM file
#####1. 处理.sam文件,提取mutation的tag,并转成.bed格式
perl /Share/home/zhangqf7/gongjing/software/ctk-master/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file iCLIP.mutation.txt ../iCLIP.sam iCLIP.tag.bed # uniq + multi
# perl /Share/home/zhangqf7/gongjing/software/ctk-master/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file iCLIP.mutation.txt ../uniqmap.sam iCLIP.tag.bed # uniq only


####2. Collapsing PCR duplicates
####2. 根据Mapping的位置以及reads的barcode去重复
perl /Share/home/zhangqf7/gongjing/software/ctk-master/tag2collapse.pl -v -big --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name iCLIP.tag.bed iCLIP.tag.uniq.bed 
perl /Share/home/zhangqf7/gongjing/software/ctk-master/selectRow.pl -q 3 -f 3 iCLIP.mutation.txt iCLIP.tag.uniq.bed > iCLIP.tag.uniq.mutation.txt 


###3. After getting the unique tags of each library, one might concatenate biological replicates, which are distinguished by different colors
###3. 作图用,将一个样本的tag用一种颜色标记
perl /Share/home/zhangqf7/gongjing/software/ctk-master/bed2rgb.pl -v -col "128,0,0" iCLIP.tag.uniq.bed iCLIP.tag.uniq.rgb.bed 


###4. As a diagnostic step, get the length distribution of unique tags, which should be a more faithful representation of the library:
###4. 统计各种长度的tag分布频数,以评估测序片段的长度
awk '{print $3-$2}' iCLIP.tag.bed  | sort -n | uniq -c | awk '{print $2"\t"$1}' > iCLIP.uniq.len.dist.txt

###由于zebrafish的注释很多类似于Zv9_NA101的name，容易造成后期分析时tags数与mutation的name不匹配，因此直接去掉
### 在我们这里是去掉chrUn_的染色体
# grep -v 'Zv9' iCLIP.tag.bed >iCLIP.tag.bed2
# grep -v 'Zv9' iCLIP.mutation.txt >iCLIP.mutation.txt2
grep -vE "chrUn|chrM" iCLIP.tag.bed > iCLIP.tag.uniq.nochrUn.bed
grep -vE "chrUn|chrM" iCLIP.mutation.txt > iCLIP.tag.uniq.mutation.nochrUn.txt
# awk '{if($9==">") {print $0}}' iCLIP.mutation.nochrUn.txt | cut -f 1-6 > iCLIP.tag.uniq.sub.bed
# awk '{if($9=="-") {print $0}}' iCLIP.mutation.nochrUn.txt | cut -f 1-6 > iCLIP.tag.uniq.del.bed
# awk '{if($9=="+") {print $0}}' iCLIP.mutation.nochrUn.txt | cut -f 1-6 > iCLIP.tag.uniq.ins.bed
perl /Share/home/zhangqf7/gongjing/software/ctk-master/getMutationType.pl -t del iCLIP.tag.uniq.mutation.nochrUn.txt iCLIP.tag.uniq.del.bed
perl /Share/home/zhangqf7/gongjing/software/ctk-master/getMutationType.pl -t ins iCLIP.tag.uniq.mutation.nochrUn.txt iCLIP.tag.uniq.ins.bed
perl /Share/home/zhangqf7/gongjing/software/ctk-master/getMutationType.pl -t sub iCLIP.tag.uniq.mutation.nochrUn.txt iCLIP.tag.uniq.sub.bed

# grep -v "chrUn" iCLIP.tag.bed > iCLIP.tag.nochrUn.bed
# grep -v "chrUn" iCLIP.mutation.txt

########################################################
##CIMS（无先验知识，mutation模型和truncation模型都予以采用）
########################################################
cd /Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/${sample}/bwa/${CTK_N}
mkdir -p CIMS
cd CIMS

###Mutation Mode
perl /Share/home/zhangqf7/gongjing/software/ctk-master/CIMS.pl -n 10 -p -v --keep-cache -c cache_mut ../iCLIP.tag.uniq.nochrUn.bed ../iCLIP.tag.uniq.sub.bed iCLIP.tag.uniq.mut.CIMS.txt

awk '{if($9<=0.05) {print $0}}' iCLIP.tag.uniq.mut.CIMS.txt|sort -k 9,9n -k 8,8nr -k 7,7n > iCLIP.tag.uniq.mut.CIMS.p05.txt
cut -f 1-6 iCLIP.tag.uniq.mut.CIMS.p05.txt > iCLIP.tag.uniq.mut.CIMS.p05.bed
# awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t"$5"\t"$6}' iCLIP.tag.uniq.mut.CIMS.p05.bed > iCLIP.tag.uniq.mut.CIMS.p05.21nt.bed 



########################################################
##CITS
########################################################
cd /Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/${sample}/bwa/${CTK_N}
mkdir -p CITS
cd CITS

perl /Share/home/zhangqf7/gongjing/software/ctk-master/removeRow.pl -q 3 -f 3 -v ../iCLIP.tag.uniq.nochrUn.bed ../iCLIP.tag.uniq.del.bed > iCLIP.tag.uniq.clean.bed
perl /Share/home/zhangqf7/gongjing/software/ctk-master/bedExt.pl -n up -l "-1" -r "-1" -v iCLIP.tag.uniq.clean.bed iCLIP.tag.uniq.clean.trunc.bed  

perl /Share/home/zhangqf7/gongjing/software/ctk-master/tag2cluster.pl -big -s -maxgap "-1" -of bed -v ../iCLIP.tag.bed  iCLIP.tag.uniq.cluster.0.bed
awk '{if($5>2) {print $0}}' iCLIP.tag.uniq.cluster.0.bed > iCLIP.tag.uniq.cluster.bed 

perl /Share/home/zhangqf7/gongjing/software/ctk-master/tag2peak.pl -big -ss -v --prefix "CITS" -gap 25 -p 0.05 -gene iCLIP.tag.uniq.cluster.bed iCLIP.tag.uniq.clean.trunc.bed iCLIP.tag.uniq.clean.CITS.p05.bed 


##################################
#   READS分布情况,tag 
##################################
# cd /Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/${sample}/bwa
# awk 'NR>1134 && $5>20 {print "chr"$3"\t"$4"\t"($4+length($10))}' uniqmap.sam| intersectBed -a - -b /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Zebrafish/ENSEMBLE_79/GeneList/single-transcript-info-mRNA -wa -wb > reads.annotation
# sed 's/_/\t/g' reads.annotation | cut -f12 | sort |uniq -c | awk '{sub(/^ */,"");sub(/ *$/,"")}1' | sed 's/ /\t/g' | awk '{print $2"\t"$1}'> reads.annotation.count

# cd ${data}/bwa/CTK_Procedure
# awk '{print "chr"$1"\t"$2"\t"$3"\t"$6}' iCLIP.tag.bed | intersectBed -a - -b /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Zebrafish/ENSEMBLE_79/GeneList/single-transcript-info-mRNA -wa -wb|awk '$4==$8'|awk '!a[$0]++' > tag.annotation
# sed 's/_/\t/g' tag.annotation | cut -f13 | sort |uniq -c | awk '{sub(/^ */,"");sub(/ *$/,"")}1' | sed 's/ /\t/g' | awk '{print $2"\t"$1}'> tag.annotation.count


# bedtools annotate -counts -i iCLIP.tag.uniq.clean.CITS.p05.bed -files /Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/danRerV10.refGene.bed6.CDS /Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/danRerV10.refGene.bed6.UTR5 /Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/danRerV10.refGene.bed6.UTR3 /Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/danRerV10.refGene.bed6.intron > iCLIP.tag.uniq.clean.CITS.p05.bed.anno
