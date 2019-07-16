#! /bin/bash

for i in egg 1cell 4cell 64cell 1K sphere shield epiboly
do
    bgfile=DMSO_$i"_background"
    fgfile=NAI_$i"_forground"

    bsub -q Z-ZQF -n 10 -oo $i.pipeline.out -eo $i.pipeline.err "python /Share/home/zhangqf7/jinsong_zhang/icSHAPE/icSHAPE-master/scripts/normalize_calcenrich_pipeline.py -f $fgfile -b $bgfile -p $i"
done
