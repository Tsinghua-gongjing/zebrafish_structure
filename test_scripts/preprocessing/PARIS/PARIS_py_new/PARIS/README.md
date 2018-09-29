# PARIS Analysis Pipeline

#### Dependency:
* pysam package
* g++ (version >= 4.9)
* bedtools

#### Installation

```bash
make
python paris.py				# Run with c++ code(fast)
python paris-old.py		# Run with python code(slow)
```

#### Recommend Parameters

```bash
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
```
--intronAnnoFile is a file produced using scripts/parseBedFromGTF.py

```
scripts/parseBedFromGTF.py -g GRCh38.gtf -o GRCh38 -s gencode
```
It produced GRCh38.genomeCoor.bed and GRCh38.transCoor.bed. Then use GRCh38.genomeCoor.bed as the parameter of --intronAnnoFile 

prefix.genomeCoor.bed looks like this:

```
chr1    29554   31097   +       MIR1302-2HG_ENSG00000243485.5   ENST00000473358.1       lincRNA 29554-30039,30564-30667,30976-31097
chr1    30267   31109   +       MIR1302-2HG_ENSG00000243485.5   ENST00000469289.1       lincRNA 30267-30667,30976-31109
chr1    30366   30503   +       MIR1302-2_ENSG00000284332.1     ENST00000607096.1       miRNA   30366-30503
chr1    34554   36081   -       FAM138A_ENSG00000237613.2       ENST00000417324.1       lincRNA 35721-36081,35277-35481,34554-35174
chr1    35245   36073   -       FAM138A_ENSG00000237613.2       ENST00000461467.1       lincRNA 35721-36073,35245-35481
```

prefix.transCoor.bed looks like this: 
```
ENST00000473358.1       MIR1302-2HG_ENSG00000243485.5   lincRNA 712     1-486,487-590,591-712
ENST00000469289.1       MIR1302-2HG_ENSG00000243485.5   lincRNA 535     1-401,402-535
ENST00000607096.1       MIR1302-2_ENSG00000284332.1     miRNA   138     1-138
ENST00000417324.1       FAM138A_ENSG00000237613.2       lincRNA 1187    1-361,362-566,567-1187
ENST00000461467.1       FAM138A_ENSG00000237613.2       lincRNA 590     1-353,354-590
ENST00000606857.1       OR4G4P_ENSG00000268020.3        unprocessed_pseudogene  840     1-840
ENST00000492842.1       OR4G11P_ENSG00000240361.1       unprocessed_pseudogene  940     1-940
ENST00000335137.3       OR4F5_ENSG00000186092.4 protein_coding  918     1-918   916-918
ENST00000466430.5       RP11-34P13.7_ENSG00000238009.6  lincRNA 2748    1-158,159-263,264-413,414-2748
```


#### Example

```bash
samtools view -bh DG.sam | samtools sort > DG.bam
samtools index DG.sam
```
<img src="https://ws4.sinaimg.cn/large/006tNbRwly1fhc89cjzkdj31010rq77b.jpg">



#### Other Scripts
##### 1. DG2Frame.py
This is An Example:

```bash
ROOT=/150T/zhangqf/tmp_lp/PARIS/tanglei
python DG2Frame.py \
    -i $ROOT/59_24/59_24_DG \
    -g ~/lipan/DYNAMIC/GTF/Gencode/GRCh38.gtf \
    -o $ROOT/59_24/analysis/59_24 \
    -s gencode \
    --genomeCoor ~/lipan/paris/human_virus/GRCh38.genomeCoor.bed
```
Produced 3 files: 

* `$ROOT/59_24/analysis/59_24.Frame`

```
#Group  lchr    lstrand lstart  lend    rchr    rstrand rstart  rend    support lcount  rcount  score
7747163 KU501215.1      +       3408    3425    chr1    -       148038821       148038844       2       6177    15857   0.000303
7747165 KU501215.1      +       3408    3436    chr1    -       172993707       172993722       2       6177    2       0.022036
7747173 KU501215.1      +       3411    3428    chr1    -       34377046        34377063        2       6177    2       0.022036
```

* `$ROOT/59_24/analysis/59_24.gFrame`

```
# From: /150T/zhangqf/tmp_lp/PARIS/tanglei/59_24/analysis/59_24.frame; Date: 2017-07-12 20:33:56.876354; Annotation: /Share/home/zhangqf8/lipan/DYNAMIC/GTF/Gencode/GRCh38.gtf
#Group  lchr    lstrand lstart  lend    rchr    rstrand rstart  rend    support lcount  rcount  score   lgenes  rgenes
7747163 KU501215.1      +       3408    3425    chr1    -       148038821       148038844       2       6177    15857   0.000303        NULL    chr1:148038822-148038844:ENSG00000206585.1:RNVU1-7:73-95
7747165 KU501215.1      +       3408    3436    chr1    -       172993707       172993722       2       6177    2       0.022036        NULL    NULL
7747173 KU501215.1      +       3411    3428    chr1    -       34377046        34377063        2       6177    2       0.022036        NULL    NULL
```

* `$ROOT/59_24/analysis/59_24.tFrame`

```
# From: /150T/zhangqf/tmp_lp/PARIS/tanglei/59_24/analysis/59_24.frame; Date: 2017-07-12 20:35:19.706766; Annotation: /Share/home/zhangqf8/lipan/DYNAMIC/GTF/Gencode/GRCh38.gtf
#Group  lchr    lstrand lstart  lend    rchr    rstrand rstart  rend    support lcount  rcount  score   ltrans  rtrans
7747163 KU501215.1      +       3408    3425    chr1    -       148038821       148038844       2       6177    15857   0.000303        NULL    chr1:148038822-148038844:ENST00000383858.1:73-95:RNVU1-7:snRNA:1-164:164
7747165 KU501215.1      +       3408    3436    chr1    -       172993707       172993722       2       6177    2       0.022036        NULL    NULL
7747173 KU501215.1      +       3411    3428    chr1    -       34377046        34377063        2       6177    2       0.022036        NULL    NULL
```



#### Bug to Fix

1. PrintSopportSam sequence and quanlity format
