#! /bin/bash

##functions##
function job_id()
{
    echo $2 | gawk 'match($0, /<([0-9]+)>/, id){print id[1]}'
}

##PIPELINE PARAMETERS###

##USER defined input file and working dir
LIBRARY_NAME=$1  ##assume that the input file in a single [LIBRARY_NAME].fastq file
SOURCEDIR=$2 ## where is the [LIBRARY_NAME].fastq
WD=$3 ##where to write the output

##current fix parameter
ADAPTOR_3="AGATCGGAAGAGCACACGTCTG" #3'adaptor sequence
MINLENGTH=25 #min length of read to report
UMILEN=13 #leading umi length for headcrop
TRANSCRIPTOME_INDEX=/Share/home/zhangqf7/jinsong_zhang/genome/zebrafish/z10/transcriptome_refSeq/danRer10.refSeq.transcriptome 
MINRPKM=1 #filter out low expression transcripts
##arg check
if [ $# != 3 ]; then
        echo "Usage: icshape_pipeline library_name sourcedir outputdir";
        exit;
    fi
##simple check
if [ ! -f "$SOURCEDIR/$LIBRARY_NAME.fastq" ];then
   echo "$SOURCEDIR/$LIBRARY_NAME.fastq dosen't exist!"
   exit
fi
##ICSHAPE SCRIPT DIR##
BIN="/Share/home/zhangqf7/jinsong_zhang/icSHAPE/icSHAPE-master/scripts"
#source /Share/home/zhangqf7/jinsong_zhang/.bash_profile_zjs 

##preprocessing
#STEP0: trim adaptor 
PROCESS_DIR=$WD/preprocessing 
LOG_DIR=$WD/log
FASTQC_DIR=$WD/fastqc_result
if [ ! -d $PROCESS_DIR ]
then
   mkdir -p $PROCESS_DIR
fi
mkdir -p $LOG_DIR
mkdir -p $FASTQC_DIR

BSUB_0=$(bsub -q Z-ZQF  -n 2 -J 'ZJS_cutadaptor_'$LIBRARY_NAME -e $LOG_DIR/$LIBRARY_NAME'_cutadapt.err' -o $LOG_DIR/$LIBRARY_NAME'_cutadapt.out' \
  "cutadapt -a $ADAPTOR_3 -m $MINLENGTH $SOURCEDIR/$LIBRARY_NAME.fastq -o $PROCESS_DIR/$LIBRARY_NAME.trimmed.fastq --untrimmed-output $PROCESS_DIR/$LIBRARY_NAME.untrimmed.fastq 1>$PROCESS_DIR/$LIBRARY_NAME.cutadapt.stdout 2>$PROCESS_DIR/$LIBRARY_NAME.cutadapt.stderr")

#STEP1: quality trim

BSUB_1=$(bsub -q Z-ZQF  -n 20 -J 'ZJS_qualtrim_'$LIBRARY_NAME -w "done($(job_id $BSUB_0))" -e $LOG_DIR/$LIBRARY_NAME'_qualtrim.err' -o $LOG_DIR/$LIBRARY_NAME'_qualtrim.out' \
  "java -Xmx30g -jar /Share/home/zhangqf/shaodi/app/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads 20 -phred33 $PROCESS_DIR/$LIBRARY_NAME.trimmed.fastq $PROCESS_DIR/$LIBRARY_NAME.trimmed.qual.fastq TRAILING:20 MINLEN:$MINLENGTH 2> $PROCESS_DIR/$LIBRARY_NAME.TRIMMOMATIC.err")

bsub -q Z-ZQF  -n 5 -J 'ZJS_fastqc_trim_'$LIBRARY_NAME -w "done($(job_id $BSUB_1))" -e $LOG_DIR/$LIBRARY_NAME'_fastqc_trim.err' -o $LOG_DIR/$LIBRARY_NAME'_fastqc_trim.out' \
  "fastqc -o $FASTQC_DIR -t 5 $PROCESS_DIR/$LIBRARY_NAME.trimmed.qual.fastq"
#STEP2: collapse duplicates
BSUB_2=$(bsub -q Z-ZQF  -n 10 -J 'ZJS_collapsing_'$LIBRARY_NAME -w "done($(job_id $BSUB_1))" -e $LOG_DIR/$LIBRARY_NAME'_collapsing.err' -o $LOG_DIR/$LIBRARY_NAME'_collapsing.out' \
  "$BIN/readCollapse.pl -U $PROCESS_DIR/$LIBRARY_NAME.trimmed.qual.fastq -o $PROCESS_DIR/$LIBRARY_NAME.rmdup.trimmed.qual.fastq 1>$PROCESS_DIR/$LIBRARY_NAME.collapsing.stdout 2>$PROCESS_DIR/$LIBRARY_NAME.collapsing.stderr")

bsub -q Z-ZQF  -n 5 -J 'ZJS_fastqc_collapsing_'$LIBRARY_NAME -w "done($(job_id $BSUB_2))" -e $LOG_DIR/$LIBRARY_NAME'_fastqc_collapse.err' -o $LOG_DIR/$LIBRARY_NAME'_fastqc_collapse.out' \
  "fastqc -o $FASTQC_DIR -t 5 $PROCESS_DIR/$LIBRARY_NAME.rmdup.trimmed.qual.fastq"

#STEP3: cut leading molecular barcode

BSUB_3=$(bsub -q Z-ZQF  -n 2 -J 'ZJS_headcropping_'$LIBRARY_NAME -w "done($(job_id $BSUB_2))" -e $LOG_DIR/$LIBRARY_NAME'_headcropping.err' -o $LOG_DIR/$LIBRARY_NAME'_headcropping.out' \
  "cutadapt -u $UMILEN -m $MINLENGTH $PROCESS_DIR/$LIBRARY_NAME.rmdup.trimmed.qual.fastq -o $PROCESS_DIR/$LIBRARY_NAME.clean.fastq 1>$PROCESS_DIR/$LIBRARY_NAME.headcropping.stdout 2>$PROCESS_DIR/$LIBRARY_NAME.headcropping.stderr")

bsub -q Z-ZQF -n 5 -J 'ZJS_fastqc_clean_'$LIBRARY_NAME -w "done($(job_id $BSUB_3))" -e $LOG_DIR/$LIBRARY_NAME'_fastqc_clean.err' -o $LOG_DIR/$LIBRARY_NAME'_fastqc_clean.out' \
  "fastqc -o $FASTQC_DIR -t 5 $PROCESS_DIR/$LIBRARY_NAME.clean.fastq"
#STEP4: mapping
MAPPING_DIR=$WD/mapping
if [ ! -d $MAPPING_DIR ]
then
    mkdir -p $MAPPING_DIR
fi
#mkdir -p $MAPPING_DIR/ribosomalRNA
mkdir -p $MAPPING_DIR/transcriptome

##mapping to transcriptome
BSUB_4=$(bsub -q Z-ZQF  -n 20 -J 'ZJS_mapping_'$LIBRARY_NAME -w "done($(job_id $BSUB_3))" -e $LOG_DIR/$LIBRARY_NAME'_mapping_transcriptome.err' -o $LOG_DIR/$LIBRARY_NAME'_mapping_transcriptome.out' \
 "bowtie2 -U $PROCESS_DIR/$LIBRARY_NAME.clean.fastq -S $MAPPING_DIR/transcriptome/$LIBRARY_NAME.transcriptome.sam -x $TRANSCRIPTOME_INDEX --no-unal --non-deterministic --time -p 20 \
  1>$MAPPING_DIR/transcriptome/$LIBRARY_NAME.transcriptome.stdout 2>$MAPPING_DIR/transcriptome/$LIBRARY_NAME.transcriptome.log")

bsub -q Z-ZQF  -J 'ZJS_rm_'$LIBRARY_NAME -w "done($(job_id $BSUB_2))" "rm -r $PROCESS_DIR/tmp"
##generate results end at rt file
RESULT_DIR=$WD/results/
if [ ! -d $RESULT_DIR ]
then
   mkdir -p $RESULT_DIR
fi
mkdir -p $RESULT_DIR/rpkm
mkdir -p $RESULT_DIR/rt
##estimateRPKM
BSUB_5=$(bsub -q Z-ZQF  -n 5 -J 'ZJS_estimateRPKM_'$LIBRARY_NAME -w "done($(job_id $BSUB_4))" -e $LOG_DIR/$LIBRARY_NAME'_estRPKM.err' -o $LOG_DIR/$LIBRARY_NAME'_estRPKM.out' \
  "$BIN/estimateRPKM.pl -i $MAPPING_DIR/transcriptome/$LIBRARY_NAME.transcriptome.sam -o $RESULT_DIR/rpkm/$LIBRARY_NAME.transcriptome.rpkm")

##calcRT
BSUB_6=$(bsub -q Z-ZQF  -n 10 -J 'ZJS_calcRT_'$LIBRARY_NAME -w "done($(job_id $BSUB_5))" -e $LOG_DIR/$LIBRARY_NAME'_calcRT.err' -o $LOG_DIR/$LIBRARY_NAME'_calcRT.out' \
 "$BIN/calcRT.pl -i $MAPPING_DIR/transcriptome/$LIBRARY_NAME.transcriptome.sam -o $RESULT_DIR/rt/$LIBRARY_NAME.transcriptome.rt -r $RESULT_DIR/rpkm/$LIBRARY_NAME.transcriptome.rpkm -c $MINRPKM")
#rm -r $PROCESS_DIR/tmp
