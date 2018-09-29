


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   Recommended Parameters
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

ANNOTATION=/Share/home/zhangqf8/lipan/paris/tanglei_test/annotation
DATA=/Share/home/zhangqf8/lipan/paris/Cell2017_data/SRR2814763/test_small_sample
OUT=$DATA

bsub -q Z-ZQF \
    ~/lipan/paris/Cell2017_data/SRR2814763/test_small_sample/paris_c++/paris-new.py \
    -i $DATA/Aligned.out.sam \
    -j $DATA/Chimeric.out.junction \
    -s $DATA/Chimeric.out.sam \
    -o $OUT/DG \
    --genomeFile $ANNOTATION/human-virus.fa \
    --intronAnnoFile ~/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed \
    --tmpFileName 2017_7_7 \
    --log $OUT/log.2017_7_7.txt \
    --error $OUT/error.2017_7_7.txt \
    --removeRedundancy yes \
    --minOverhang 15 \
    --localAlign no \
    --preserveMultimap no \
    --intronFlanking 3 \
    --minOverlap 5 \
    --multipleDG no \
    --maxGap 10 \
    --maxDGOverhang 30 \
    --coverage pileup \
    --genomeSizeFile $ANNOTATION/genome.size \
    --minSupport 2 \
    --scoringMethod harmonic







# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   Other Test Scripts
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




python script/main.py \
    -i /150T/zhangqf/tmp_lp/PARIS/766_72/Aligned.out.sam \
    -j /150T/zhangqf/tmp_lp/PARIS/766_72/Chimeric.out.junction \
    -s /150T/zhangqf/tmp_lp/PARIS/766_72/Chimeric.out.sam \
    -o /Share/home/zhangqf8/lipan/paris/test_python/DG \
    -a /Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed \
    -g /Share/home/zhangqf8/lipan/paris/debug_test/annotation/human-virus.fa \
    -z /Share/home/zhangqf8/lipan/paris/debug_test/annotation/genome.size \
    --name big


# test in /Share/home/zhangqf8/lipan/paris/test_python
SCRIPTS=/Share/home/zhangqf8/lipan/paris/test_python/script
PYTHONPATH=$SCRIPTS:$PYTHONPATH
python $SCRIPTS/main.py \
    -i /Share/home/zhangqf8/lipan/paris/debug_test/rawData/Aligned.out.sort.sam \
    -j /Share/home/zhangqf8/lipan/paris/debug_test/rawData/Chimeric.out.junction \
    -s /Share/home/zhangqf8/lipan/paris/debug_test/rawData/Chimeric.out.sam \
    -o /Share/home/zhangqf8/lipan/paris/test_python/DG \
    -a /Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed \
    -g /Share/home/zhangqf8/lipan/paris/debug_test/annotation/human-virus.fa \
    -z /Share/home/zhangqf8/lipan/paris/debug_test/annotation/genome.size \
    --name try



# test in /Share/home/zhangqf8/lipan/paris/766_72/python
SCRIPTS=/Share/home/zhangqf8/lipan/paris/python_scripts
ANNOTATION=/Share/home/zhangqf8/lipan/paris/tanglei_test/annotation
DATA=/150T/zhangqf/tmp_lp/PARIS/766_72
PYTHONPATH=$SCRIPTS:$PYTHONPATH
python $SCRIPTS/main.py \
    -i $DATA/Aligned.out.sam \
    -j $DATA/Chimeric.out.junction \
    -s $DATA/Chimeric.out.sam \
    -o /Share/home/zhangqf8/lipan/paris/766_72/python/766_72_DG \
    -a /Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed \
    -g $ANNOTATION/human-virus.fa \
    -z $ANNOTATION/genome.size \
    --name 2017_6_25 > log.2017_6_25.txt 2> error.2017_6_25.txt


sys.argv = ['bin.py', '-i', '/150T/zhangqf/tmp_lp/PARIS/766_72/Aligned.out.sam', 
    '-j', '/150T/zhangqf/tmp_lp/PARIS/766_72/Chimeric.out.junction',
    '-s', '/150T/zhangqf/tmp_lp/PARIS/766_72/Chimeric.out.sam',
    '-o', '/Share/home/zhangqf8/lipan/paris/766_72/python/766_72_DG',
    '-a', '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed',
    '-g', '/Share/home/zhangqf8/lipan/paris/tanglei_test/annotation/human-virus.fa', 
    '-z', '/Share/home/zhangqf8/lipan/paris/tanglei_test/annotation/genome.size',
    '--name', '2017_6_25']


# test in /Share/home/zhangqf8/lipan/splash/fastq
SCRIPTS=~/lipan/paris/python_scripts
ANNOTATION=/Share/home/zhangqf8/lipan/DYNAMIC/GTF
DATA=/Share/home/zhangqf8/lipan/splash/fastq
export PYTHONPATH=$SCRIPTS:$PYTHONPATH
stdbuf -oL \
    python $SCRIPTS/main.py \
    -i $DATA/SRR3404926_star.Aligned.out.sam \
    -j $DATA/SRR3404926_star.Chimeric.out.junction \
    -s $DATA/SRR3404926_star.Chimeric.out.sam \
    -o $DATA/SRR3404926_paris \
    -a /Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed \
    -g $ANNOTATION/hg38.fa \
    -z $ANNOTATION/size/hg38.genome.size \
    --name 2017_6_25 \
    2> $DATA/error.SRR3404926_paris.2017_6_26.txt \
    > $DATA/log.SRR3404926_paris.2017_6_26.txt


ANNOTATION="/Share/home/zhangqf8/lipan/DYNAMIC/GTF"
DATA="/Share/home/zhangqf8/lipan/splash/fastq"
sys.argv = ['bin.py', '-i', DATA+'/SRR3404926_star.Aligned.out.sam', 
    '-j', DATA+'/SRR3404926_star.Chimeric.out.junction',
    '-s', DATA+'/SRR3404926_star.Chimeric.out.sam',
    '-o', DATA+'/test_python/SRR3404926_paris',
    '-a', '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed',
    '-g', ANNOTATION+'/hg38.fa', 
    '-z', ANNOTATION+'/size/hg38.genome.size',
    '--multi',
    '--log', DATA+'/test_python/log_paris.2017_6_27.txt',
    '--error', DATA+'/test_python/error_paris.2017_6_27.txt',
    '--removeRedundancy',
    '--name', '2017_6_25']


# run in /Share/home/zhangqf8/lipan/paris/Cell2017_data/SRR2814763/test_python

ANNOTATION="/Share/home/zhangqf8/lipan/DYNAMIC/GTF"
DATA="/Share/home/zhangqf8/lipan/paris/Cell2017_data/SRR2814763"
sys.argv = ['bin.py', '-i', DATA+'/Aligned.out.sam', 
    '-j', DATA+'/Chimeric.out.junction',
    '-s', DATA+'/Chimeric.out.sam',
    '-o', DATA+'/test_python/DG',
    '-a', '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed',
    '-g', ANNOTATION+'/hg38.fa', 
    '-z', ANNOTATION+'/size/hg38.genome.size',
    '--multi',
    '--log', DATA+'/test_python/log.2017_6_28.txt',
    '--error', DATA+'/test_python/error.2017_6_28.txt',
    '--removeRedundancy',
    '--name', '2017_6_28']



# new parameters

DATA=~/lipan/paris/Cell2017_data/SRR2814763
ANNOTATION=~/lipan/DYNAMIC/GTF
nohup ../paris.py \
    -i $DATA/Aligned.out.sam \
    -j $DATA/Chimeric.out.junction \
    -s $DATA/Chimeric.out.sam \
    -o ./DG \
    --genomeFile $ANNOTATION/hg38.fa \
    --intronAnnoFile ~/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed \
    --tmpFileName 2017_6_30 \
    --log ./log.2017_6_30.txt \
    --error ./error.2017_6_30.txt \
    --removeRedundancy yes \
    --minOverhang 15 \
    --localAlign no \
    --preserveMultimap yes \
    --intronFlanking 3 \
    --minOverlap 5 \
    --multipleDG yes \
    --maxGap 10 \
    --maxDGOverhang 30 \
    --coverage pileup \
    --genomeSizeFile $ANNOTATION/size/hg38.genome.size \
    --minSupport 2 \
    --scoringMethod harmonic &



# transcriptome work ?

DATA=/150T/zhangqf/tmp_lp/PARIS/testTranscriptome
nohup ../paris.py \
    -i $DATA/Aligned.out.sam \
    -j $DATA/Chimeric.out.junction \
    -s $DATA/Chimeric.out.sam \
    -o ./transcriptome_DG \
    --genomeFile $DATA/766wt.fa \
    --tmpFileName 2017_6_30 \
    --log ./log.2017_6_30.txt \
    --error ./error.2017_6_30.txt \
    --removeRedundancy yes \
    --minOverhang 15 \
    --localAlign no \
    --preserveMultimap yes \
    --intronFlanking 3 \
    --minOverlap 5 \
    --multipleDG yes \
    --maxGap 10 \
    --maxDGOverhang 30 \
    --coverage pileup \
    --genomeSizeFile $DATA/chrNameLength.txt \
    --minSupport 2 \
    --scoringMethod harmonic &




### test paris-new.py

ANNOTATION="/Share/home/zhangqf8/lipan/DYNAMIC/GTF"
DATA="/Share/home/zhangqf8/lipan/paris/Cell2017_data/SRR2814763/test_small_sample"
sys.argv = ['bin.py', '-i', DATA+'/Aligned.out.sam', 
    '-j', DATA+'/Chimeric.out.junction',
    '-s', DATA+'/Chimeric.out.sam',
    '-o', DATA+'/DG',
    '--intronAnnoFile', '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed',
    '--genomeFile', ANNOTATION+'/hg38.fa', 
    '--genomeSizeFile', ANNOTATION+'/size/hg38.genome.size',
    '--log', DATA+'/log.2017_6_28.txt',
    '--error', DATA+'/error.2017_6_28.txt',
    '--tmpFileName', '2017_7_6']






ANNOTATION=/Share/home/zhangqf8/lipan/paris/tanglei_test/annotation
DATA=/150T/zhangqf/tmp_lp/PARIS/huh7_1
OUT=/Share/home/zhangqf8/lipan/paris/huh7_1

bsub -q Z-ZQF ~/lipan/paris/Cell2017_data/SRR2814763/test_small_sample/paris_c++/paris-new.py \
    -i $DATA/Aligned_prim_N.out.gapped_header.sam \
    -j $DATA/Chimeric.out.junction \
    -s $DATA/Chimeric.out.sam \
    -o $OUT/DG \
    --genomeFile $ANNOTATION/human-virus.fa \
    --intronAnnoFile ~/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed \
    --tmpFileName 2017_7_7 \
    --log $OUT/log.2017_7_7.txt \
    --error $OUT/error.2017_7_7.txt \
    --removeRedundancy yes \
    --minOverhang 15 \
    --localAlign no \
    --preserveMultimap no \
    --intronFlanking 3 \
    --minOverlap 5 \
    --multipleDG no \
    --maxGap 10 \
    --maxDGOverhang 30 \
    --coverage pileup \
    --genomeSizeFile $ANNOTATION/genome.size \
    --minSupport 2 \
    --scoringMethod harmonic















