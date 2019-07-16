
### Usage
# bash paris_preprocessing.sh shi-zp-1-rep-combine-89

# output: (under dir: /Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-1-rep-combine-89)
# ./clean_data
# 	1-ori.fastq  
# 	2-cut.fastq  
# 	2-cut_filter.fastq  
# 	3-derep.fastq  
# 	4-trim.fastq  
# 	5-qua.fastq # clean reads for mapping
# ./27-DG # duplex files


### 

Trimmomatic=/Share/home/zhangqf7/gongjing/zebrafish/script/icSHAPE-master/bin/trimmomatic-0.30.jar
star=/Share/home/zhangqf/usr/bin/STAR
paris=/Share/home/zhangqf7/gongjing/zebrafish/script/PARIS_py_new/PARIS/paris.py # new copy from /tmp/test_paris_pipeline/PARIS
index=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/combine_transcriptome_star # merged transcriptome index, refseq+ensembl89
fasta=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/refseq_ensembl91_merge.tarns.fa
chrsize=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/combine_transcriptome_star91/chrNameLength.txt
adapter_fa=/Share/home/zhangqf7/gongjing/zebrafish/script/adapter.fa

sample=$1
echo $sample

thread=10

# define sample dir
cd /Share/home/zhangqf7/gongjing/zebrafish/data/paris/${sample}

mkdir  -p clean_data/fastqc 

########## cut adapter
cutadapt -a GATCGGAAGAGCACACGTCTG -n 5 -e 0.15 -m 30 -o clean_data/2-cut.fastq clean_data/1-ori.fastq
fastqc -t $thread clean_data/1-ori.fastq -o clean_data/fastqc

### filter by sequencing 
awk 'BEGIN{permision=0}{if(NR%4==1){tmp=$0}else{if(NR%4==2){if($0~/^....TG/){print tmp; print $0; permission=1}else{permission=0}}else{if(permission==1){print $0}} }}' clean_data/2-cut.fastq > clean_data/2-cut_filter.fastq
fastqc -t $thread clean_data/2-cut.fastq -o clean_data/fastqc

########## derep
perl $script/fastq2collapse.pl clean_data/2-cut_filter.fastq clean_data/3-derep.fastq 
fastqc -t $thread clean_data/2-cut_filter.fastq -o clean_data/fastqc

######### cut barcode
perl $script/stripBarcode.pl -format fastq -len 9 clean_data/3-derep.fastq clean_data/4-trim.fastq
fastqc -t $thread clean_data/3-derep.fastq  -o clean_data/fastqc

########## trim low quality base
java -Xmx4g -jar $Trimmomatic SE -phred33 -threads $thread clean_data/4-trim.fastq clean_data/5-qua.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30
# here major discard read too short (<30bp) 
fastqc -t $thread clean_data/4-trim.fastq  -o clean_data/fastqc
fastqc -t $thread clean_data/5-qua.fastq  -o clean_data/fastqc





################# generate index
$star --runMode genomeGenerate --genomeDir $index  --genomeFastaFiles $fasta --sjdbGTFfile $gtf  --sjdbOverhang 100 --runThreadN 12

################mapping
$star --runMode alignReads --genomeDir $index --readFilesIn clean_data/5-qua.fastq --outFileNamePrefix 6-star --outReadsUnmapped rm-unmapped8 --outFilterMultimapNmax 10 --outSAMattributes All --scoreGapNoncan -4 --scoreGapATAC -4 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --runThreadN $thread

################ calling duplex
python $paris -i 6-starAligned.out.sam -j 6-starChimeric.out.junction -s 6-starChimeric.out.sam -o DG --genomeFile $fasta --tmpFileName 20180518 --removeRedundancy yes --genomeSizeFile $chrsize 



# STAR generate index of refseq/ensembl merged index
# When you have a genome with large number of references, you need to scale --genomeChrBinNbits as log2(GenomeSize/NumberOfReferences)=log2(AverageReferenceLength), her log2(2761.9)
#STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ./combine_transcriptome_star --genomeFastaFiles ./refseq_ensembl_merge.tarns.fa --genomeChrBinNbits 11




















