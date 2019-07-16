dir=/Share/home/zhangqf7/gongjing/zebrafish/data/paris/$1
thread=20
star=/Share/home/zhangqf/usr/bin/STAR
paris=/Share/home/zhangqf7/gongjing/zebrafish/script/PARIS_py_new/PARIS/paris.py
index=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/combine_transcriptome_star91
fasta=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/refseq_ensembl91_merge.tarns.fa
chrsize=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/combine_transcriptome_star91/chrNameLength.txt
cd $dir


mkdir -p $dir/downsampling_N/downsampling_$3
#/Share/home/zhangqf/usr/jre1.8.0_73/bin/java -jar /Share/home/zhangqf/usr/picard-tools-2.1.1/picard.jar DownsampleSam I=6-starAligned.out.sam O=downsampling/6-starAligned.out.ds.sam P=$2
#/Share/home/zhangqf/usr/jre1.8.0_73/bin/java -jar /Share/home/zhangqf/usr/picard-tools-2.1.1/picard.jar DownsampleSam I=6-starChimeric.out.sam O=downsampling/6-starChimeric.out.ds.sam P=$2
#cp 6-starChimeric.out.junction downsampling/

#python /Share/home/zhangqf7/gongjing/zebrafish/script/paris_sam_downsampling.py 6-starAligned.out.sam 6-starChimeric.out.sam 6-starChimeric.out.junction $dir/downsampling_N/downsampling_${3} $2 $3
cd $dir/downsampling_N/downsampling_$3

#grep -v ^@ 6-starAligned.out.ds.sam| awk '{print "@"$1"\n"$10"\n+\n"$11}' > ds.fastq
#grep -v ^@ 6-starChimeric.out.ds.sam| awk '{print "@"$1"\n"$10"\n+\n"$11}' >> ds.fastq
#mkdir -p mapping
#cd mapping
#$star --runMode alignReads --genomeDir $index --readFilesIn ../ds.fastq --outFileNamePrefix 6-star --outReadsUnmapped rm-unmapped8 --outFilterMultimapNmax 10 --outSAMattributes All --scoreGapNoncan -4 --scoreGapATAC -4 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --runThreadN $thread
#mkdir -p duplex
#cd duplex

python $paris -i 6-starAligned.out.ds.sam -j 6-starChimeric.out.ds.junction -s 6-starChimeric.out.ds.sam -o 27-DG --genomeFile $fasta --tmpFileName 20180409 --removeRedundancy yes --genomeSizeFile $chrsize

