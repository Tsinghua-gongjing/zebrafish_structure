d=/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/new_mergepeaks_d10/separate
UTRCDS_bed=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/wholeTransPositionDetailUTRCDS.bed
tranSize=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.size
tranFa=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.fa
meme=/Share/home/zhangqf7/gongjing/zebrafish/script/zhangting/paris_RBP/motif_CISBP_RNA_narrow/Collapsed.used.meme

cd $d

#cat way1.bed way2.bed way3.bed way4.bed way5.bed > way12345.bed
#cat way2.bed way3.bed way4.bed way5.bed > way2345.bed

#for i in way12345.bed
for i in way5.bed
do

sort_bed=`echo $i|sed 's/bed/sorted.bed/'`
merge_sort_bed=`echo $sort_bed|sed 's/bed/merge.bed/'`
sort -k1,1 -k2,2n $i > $sort_bed
bedtools merge -i $sort_bed > $merge_sort_bed

merge_sort_anno_bed=`echo $merge_sort_bed|sed 's/bed/anno.bed/'`
intersectBed -a $merge_sort_bed -b $UTRCDS_bed -wa -wb -f 0.5 > $merge_sort_anno_bed

merge_sort_anno_utr3_bed=`echo $merge_sort_anno_bed|sed 's/bed/utr3.bed/'`
awk '{if($8=="utr3")print $1"\t"$2"\t"$3}' $merge_sort_anno_bed > $merge_sort_anno_utr3_bed

merge_sort_anno_utr3_fa=`echo $merge_sort_anno_utr3_bed|sed 's/bed/fa/'`
bedtools getfasta -fi $tranFa -bed $merge_sort_anno_utr3_bed -fo $merge_sort_anno_utr3_fa

shuffleBed -seed 1234 -i $merge_sort_anno_utr3_bed -chrom -g $tranSize > bg/$merge_sort_anno_utr3_bed
bedtools getfasta -fi $tranFa -bed bg/$merge_sort_anno_utr3_bed -fo bg/$merge_sort_anno_utr3_fa

dreme=${merge_sort_anno_utr3_fa}_dreme
/Share/home/zhangqf7/bin/dreme -p $merge_sort_anno_utr3_fa -n bg/$merge_sort_anno_utr3_fa -norc -mink 5 -rna -o $dreme

/Share/home/zhangqf7/bin/tomtom -xalph -no-ssc -oc ./tomtom_5_Collapsed-used -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 5 -norc ./dreme.txt $meme 

done
