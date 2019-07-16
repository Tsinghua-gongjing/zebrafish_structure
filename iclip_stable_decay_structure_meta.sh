bed1=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-5/bwa/CTK_Procedure/CITS/iCLIP.tag.uniq.clean.CITS.p05.ext20.tranIn.len41.sort.merge.anno.utr3.fimo/fimo.new.HUR_motif2.stable.bed
bed2=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-5/bwa/CTK_Procedure/CITS/UTR3_minus_iCLIP.fimo/fimo.new.HUR_motif2.stable.bed
shape_out=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_iCLIP_notiCLIP_motif_structure_stable.ext20search.3.pdf

# python bed_motif_structure_meta.py $bed1:$bed2 in_iCLIP_stable:not_in_iCLIP_stable $shape_out:$shape_out $savefn 20


bed1=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-6/bwa/CTK_Procedure/CITS/iCLIP.tag.uniq.clean.CITS.p05.ext20.tranIn.len41.sort.merge.anno.utr3.fimo/fimo.new.HUR_motif2.stable.bed
bed2=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-6/bwa/CTK_Procedure/CITS/UTR3_minus_iCLIP.fimo/fimo.new.HUR_motif2.stable.bed
shape_out=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_iCLIP_notiCLIP_motif_structure_stable.ext20search.3.pdf

# python bed_motif_structure_meta.py $bed1:$bed2 in_iCLIP_stable:not_in_iCLIP_stable $shape_out:$shape_out $savefn 20

# h4 bind h6 not bind region
h4_bed=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-5/bwa/CTK_Procedure/CITS/iCLIP.tag.uniq.clean.CITS.p05.ext20.tranIn.len41.sort.merge.anno.utr3.fimo/fimo.new.HUR_motif2.bed
h6_bed=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-6/bwa/CTK_Procedure/CITS/iCLIP.tag.uniq.clean.CITS.p05.ext20.tranIn.len41.sort.merge.anno.utr3.fimo/fimo.new.HUR_motif2.bed
h4_minus_h6_bed=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.bed
# bedtools subtract -a $h4_bed -b $h6_bed|awk '{if($3-$2==7)print}' > $h4_minus_h6_bed

h6_minus_h4_bed=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_bed_minus_h4_bed.bed
 # bedtools subtract -a $h6_bed -b $h4_bed|awk '{if($3-$2==7)print}' > $h6_minus_h4_bed

# h4 bind h6 not bind region, no h6 bind tx
# awk 'NR==FNR{a[$1]=$0;next} !($1 in a){print $0}' $h6_bed h4_bed_minus_h6_bed.bed > h4_bed_minus_h6_bed_minus_H6_tx.bed

# decay gene list
# cd /Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_paper/DEGseq
# awk -F "\t" '{if($30>=0)print $2}' Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.txt|grep -v "Transcript" >Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.txt
# awk -F "\t" '{if($6>=0)print $2}' H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.txt|grep -v "Transcript" > H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.txt

# awk -F "\t" '{if($30>=0.585)print $2}' Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.txt|grep -v "Transcript" >Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.585.txt
# awk -F "\t" '{if(($30<0.585)&&($30>-0.585))print $2}' Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.txt|grep -v "Transcript" >Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.abslt0.585.txt

# decay and rescue gene list
# awk -F "\t" '{if(($30>=0)&&($31>=0))print $2}' Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.txt|grep -v "Transcript" >Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.txt
# awk -F "\t" '{if(($6>=0)&&($15>=0))print $2}' H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.txt|grep -v "Transcript" > H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0_yge0.txt

# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' /Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_paper/DEGseq/Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.txt h4_bed_minus_h6_bed.bed  > h4_bed_minus_h6_bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.bed
# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' /Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_paper/DEGseq/H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.txt h4_bed_minus_h6_bed.bed  > h4_bed_minus_h6_bed.H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.bed

# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' /Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_paper/DEGseq/Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.txt h4_bed_minus_h6_bed.bed > h4_bed_minus_h6_bed.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.bed

# h4 bind h6 not bind region, decay gene structure
bed1=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.bed
bed2=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.bed
shape_out1=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out
shape_out2=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.bed.structure.pdf
# python bed_motif_structure_meta.py $bed1:$bed2 H4:H6 $shape_out1:$shape_out2 $savefn 20

bed1=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.bed
bed2=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.bed
shape_out1=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out
shape_out2=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.bed.structure.pdf
# python bed_motif_structure_meta.py $bed1:$bed2 H4:H6 $shape_out1:$shape_out2 $savefn 20

bed1=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.585.bed
bed2=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.585.bed
shape_out1=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out
shape_out2=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.585.bed.structure.pdf
# python bed_motif_structure_meta.py $bed1:$bed2 H4:H6 $shape_out1:$shape_out2 $savefn 20

bed1=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.abslt0.585.bed
bed2=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.abslt0.585.bed
shape_out1=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out
shape_out2=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.abslt0.585.bed.structure.pdf
# python bed_motif_structure_meta.py $bed1:$bed2 H4:H6 $shape_out1:$shape_out2 $savefn 20

# h4 bind h6 not bind region, no h6 bind tx, decay gene structure

# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' /Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_paper/DEGseq/Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.txt h4_bed_minus_h6_bed_minus_H6_tx.bed  > h4_bed_minus_h6_bed_minus_H6_tx.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.bed
# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' /Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_paper/DEGseq/H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.txt h4_bed_minus_h6_bed_minus_H6_tx.bed  > h4_bed_minus_h6_bed_minus_H6_tx.bed.H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.bed

bed1=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed_minus_H6_tx.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.bed
bed2=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed_minus_H6_tx.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.bed
shape_out1=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out
shape_out2=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed_minus_H6_tx.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.bed.structure.pdf
# python bed_motif_structure_meta.py $bed1:$bed2 H4:H6 $shape_out1:$shape_out2 $savefn 20

bed1=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed_minus_H6_tx.bed.H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.bed
bed2=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed_minus_H6_tx.bed.H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.bed
shape_out1=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out
shape_out2=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed_minus_H6_tx.bed.H6_Elavl1a-MO_Control-MO_vs_H6_Elavl1a-MO_2Elavl1a-Rescue.ge0.bed.structure.pdf
# python bed_motif_structure_meta.py $bed1:$bed2 H4:H6 $shape_out1:$shape_out2 $savefn 20

# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' /Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_New_refseq/stable.txt h4_bed_minus_h6_bed.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.bed > h4_bed_minus_h6_bed.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.stable.bed

bed1=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.stable.bed
bed2=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.stable.bed
shape_out1=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out
shape_out2=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.stable.bed.structure.pdf
# python bed_motif_structure_meta.py $bed1:$bed2 H4:H6 $shape_out1:$shape_out2 $savefn 20

bed3=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.decay.bed
bed4=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.decay.bed
shape_out3=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out
shape_out4=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.decay.bed.structure.pdf
# python bed_motif_structure_meta.py $bed3:$bed4 H4:H6 $shape_out3:$shape_out4 $savefn 20

savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.bed.Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.stable_decay.bed.structure.pdf
# python bed_motif_structure_meta.py $bed1:$bed2:$bed3:$bed4 H4:H6:H4:H6 $shape_out1:$shape_out2:$shape_out3:$shape_out4 $savefn 20


shape_sphere=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out
shape_shield=/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.structure.e0.txt
# python bed_structure_compare.py $h4_minus_h6_bed H4:H6 $shape_sphere:$shape_shield $savefn 0
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.structure.e3.txt
# python bed_structure_compare.py $h4_minus_h6_bed H4:H6 $shape_sphere:$shape_shield $savefn 3
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.structure.e10.txt
# python bed_structure_compare.py $h4_minus_h6_bed H4:H6 $shape_sphere:$shape_shield $savefn 10

savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_bed_minus_h4_bed.structure.e0.txt
# python bed_structure_compare.py $h6_minus_h4_bed H4:H6 $shape_sphere:$shape_shield $savefn 0
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_bed_minus_h4_bed.structure.e3.txt
# python bed_structure_compare.py $h6_minus_h4_bed H4:H6 $shape_sphere:$shape_shield $savefn 3
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_bed_minus_h4_bed.structure.e10.txt
# python bed_structure_compare.py $h6_minus_h4_bed H4:H6 $shape_sphere:$shape_shield $savefn 10

morestructure=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.structure.e0.txt.morestructure
morestructure2=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.structure.e0.txt.morestructure2
down_regulated=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_paper/DEGseq/Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.txt
decay=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_New_refseq/decay_id.txt
savefn1=${morestructure}.KO.decay.png
savefn2=${morestructure2}.KO.decay.png
# python venn_plot.py $morestructure:$down_regulated:$decay iCLIP_more_structure:down_regulated:decay $savefn1
# python venn_plot.py $morestructure2:$down_regulated:$decay iCLIP_more_structure:down_regulated:decay $savefn2

lessstructure=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_bed_minus_h4_bed.structure.e0.txt.lessstructure
lessstructure2=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_bed_minus_h4_bed.structure.e0.txt.lessstructure2
stable=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_New_refseq/stable_id.txt
savefn3=${lessstructure}.KO.stable.png
savefn4=${lessstructure2}.KO.stable.png
# python venn_plot.py $lessstructure:$down_regulated:$stable iCLIP_less_structure:down_regulated:stable $savefn3
# python venn_plot.py $lessstructure2:$down_regulated:$stable iCLIP_less_structure:down_regulated:stable $savefn4

savefn_all=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_bed_minus_h6_bed.structure.e0.overlap.pdf
# montage -mode concatenate -tile 2x $savefn1 $savefn2 $savefn3 $savefn4 $savefn_all

h4_iclip_bed_all=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-5/bwa/CTK_Procedure/CITS/iCLIP.tag.uniq.clean.CITS.p05.ext20.tranIn.len41.sort.merge.anno.utr3.bed
h6_iclip_bed_all=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-6/bwa/CTK_Procedure/CITS/iCLIP.tag.uniq.clean.CITS.p05.ext20.tranIn.len41.sort.merge.anno.utr3.bed
h4_all_minus_h6_all=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_all_minus_h6_all
# bedtools subtract -a $h4_iclip_bed_all -b $h6_iclip_bed_all|awk '{if($3-$2>=10)print}' > $h4_all_minus_h6_all
# python bed_structure_compare.py $h4_all_minus_h6_all H4:H6 $shape_sphere:$shape_shield ${h4_all_minus_h6_all}.structure.e0.txt 0

h6_all_minus_h4_all=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_all_minus_h4_all
# bedtools subtract -a $h6_iclip_bed_all -b $h4_iclip_bed_all|awk '{if($3-$2>=10)print}' > $h6_all_minus_h4_all
# python bed_structure_compare.py $h6_all_minus_h4_all H4:H6 $shape_sphere:$shape_shield ${h6_all_minus_h4_all}.structure.e0.txt 0

morestructure=${h4_all_minus_h6_all}.structure.e0.txt.morestructure
morestructure2=${h4_all_minus_h6_all}.structure.e0.txt.morestructure2
# down_regulated=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_paper/DEGseq/Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0_yge0.txt # consider rescue
down_regulated=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_paper/DEGseq/Control_MO_vs_Elavl1_MO_vs_2Elavl1_MO_Rescue.ge0.585.txt # not consider rescue
decay=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_New_refseq/decay_id.txt
savefn1=${morestructure}.KO.decay.png
savefn2=${morestructure2}.KO.decay.png
# python venn_plot.py $morestructure:$down_regulated:$decay iCLIP_more_structure:down_regulated:decay $savefn1
# python venn_plot.py $morestructure2:$down_regulated:$decay iCLIP_more_structure:down_regulated:decay $savefn2

lessstructure=${h6_all_minus_h4_all}.structure.e0.txt.lessstructure
lessstructure2=${h6_all_minus_h4_all}.structure.e0.txt.lessstructure2
stable=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_New_refseq/stable_id.txt
zygotic=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_New_refseq/zygotic_id.txt
savefn3=${lessstructure}.KO.stable.png
savefn4=${lessstructure2}.KO.stable.png
# python venn_plot.py $lessstructure:$down_regulated:$stable iCLIP_less_structure:down_regulated:stable $savefn3
# python venn_plot.py $lessstructure2:$down_regulated:$stable iCLIP_less_structure:down_regulated:stable $savefn4

# h4_all=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_all.utr3.bed
h4_all=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_all_rep12.c6.utr3.bed
# python bed_structure_compare.py $h4_all H4:H6 $shape_sphere:$shape_shield ${h4_all}.structure.e0.txt 0
# python venn_plot.py ${h4_all}.structure.e0.txt.morestructure:$down_regulated:$decay iCLIP_more_structure:down_regulated:decay ${h4_all}.structure.e0.txt.morestructure.KO.decay.png
# python venn_plot.py ${h4_all}.structure.e0.txt.morestructure2:$down_regulated:$decay iCLIP_more_structure:down_regulated:decay ${h4_all}.structure.e0.txt.morestructure2.KO.decay.png
# python venn_plot.py ${h4_all}.structure.e0.txt.morestructure3:$down_regulated:$decay iCLIP_more_structure:down_regulated:decay ${h4_all}.structure.e0.txt.morestructure3.KO.decay.png
# python venn_plot.py ${h4_all}.structure.e0.txt.morestructure4:$down_regulated:$decay iCLIP_more_structure:down_regulated:decay ${h4_all}.structure.e0.txt.morestructure4.KO.decay.png
# python venn_plot.py ${h4_all}.structure.e0.txt.lessstructure:$down_regulated:$stable iCLIP_less_structure:down_regulated:stable ${h4_all}.structure.e0.txt.lessstructure.KO.stable.png
# python venn_plot.py ${h4_all}.structure.e0.txt.lessstructure2:$down_regulated:$stable iCLIP_less_structure:down_regulated:stable ${h4_all}.structure.e0.txt.lessstructure2.KO.stable.png
# python venn_plot.py ${h4_all}.structure.e0.txt.lessstructure3:$down_regulated:$stable iCLIP_less_structure:down_regulated:stable ${h4_all}.structure.e0.txt.lessstructure3.KO.stable.png
# python venn_plot.py ${h4_all}.structure.e0.txt.lessstructure4:$down_regulated:$stable iCLIP_less_structure:down_regulated:stable ${h4_all}.structure.e0.txt.lessstructure4.KO.stable.png
# python venn_plot.py ${h4_all}.structure.e0.txt.stablestructure:$down_regulated:$stable iCLIP_stable_structure:down_regulated:stable ${h4_all}.structure.e0.txt.stablestructure.KO.stable.png
# python venn_plot.py ${h4_all}.structure.e0.txt.stablestructure2:$down_regulated:$stable iCLIP_stable_structure:down_regulated:stable ${h4_all}.structure.e0.txt.stablestructure2.KO.stable.png
# python venn_plot.py ${h4_all}.structure.e0.txt.stablestructure3:$down_regulated:$stable iCLIP_stable_structure:down_regulated:stable ${h4_all}.structure.e0.txt.stablestructure3.KO.stable.png
# python venn_plot.py ${h4_all}.structure.e0.txt.stablestructure4:$down_regulated:$stable iCLIP_stable_structure:down_regulated:stable ${h4_all}.structure.e0.txt.stablestructure4.KO.stable.png

# python venn_plot.py ${h4_all}.structure.e0.txt.morestructure4:$decay:$stable iCLIP_more_structure:decay:stable ${h4_all}.structure.e0.txt.morestructure4.KO.decay.pdf
# python venn_plot.py ${h4_all}.structure.e0.txt.lessstructure4:$decay:$stable iCLIP_less_structure:decay:stable ${h4_all}.structure.e0.txt.lessstructure4.KO.stable.pdf
# python venn_plot.py ${h4_all}.structure.e0.txt.stablestructure4:$decay:$stable iCLIP_stable_structure:decay:stable ${h4_all}.structure.e0.txt.stablestructure4.KO.stable.pdf

# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $1}' $decay ${h4_all}.structure.e0.txt.morestructure4 > ${h4_all}.structure.e0.txt.morestructure4.decay
# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' ${h4_all}.structure.e0.txt.morestructure4.decay ${h4_all}|awk '{if($3-$2==41)print}' > ${h4_all}.structure.e0.txt.morestructure4.decay.iclip.bed 

# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $1}' $stable ${h4_all}.structure.e0.txt.lessstructure4 > ${h4_all}.structure.e0.txt.lessstructure4.stable
# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' ${h4_all}.structure.e0.txt.lessstructure4.stable ${h4_all}|awk '{if($3-$2==41)print}' > ${h4_all}.structure.e0.txt.lessstructure4.stable.iclip.bed 

# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $1}' $stable ${h4_all}.structure.e0.txt.stablestructure4 > ${h4_all}.structure.e0.txt.stablestructure4.stable
# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' ${h4_all}.structure.e0.txt.stablestructure4.stable ${h4_all}|awk '{if($3-$2==41)print}' > ${h4_all}.structure.e0.txt.stablestructure4.stable.iclip.bed 

# ol1=${h4_all}.structure.e0.txt.morestructure4.decay.iclip.bed
# ol2=${h4_all}.structure.e0.txt.lessstructure4.stable.iclip.bed
# ol3=${h4_all}.structure.e0.txt.stablestructure4.stable.iclip.bed 
# python bed_motif_structure_meta.py $ol1:$ol1 H4:H6 $shape_sphere:$shape_shield ${ol1}.structure.ext0.pdf 0 
# python bed_motif_structure_meta.py $ol2:$ol2 H4:H6 $shape_sphere:$shape_shield ${ol2}.structure.ext0.pdf 0 
# python bed_motif_structure_meta.py $ol3:$ol3 H4:H6 $shape_sphere:$shape_shield ${ol3}.structure.ext0.pdf 0 


# python venn_plot.py ${h4_all}.structure.e0.txt.morestructure4:$decay:$down_regulated iCLIP_more_structure:decay:down_regulated ${h4_all}.structure.e0.txt.morestructure4.KO.decay.down.pdf
# python venn_plot.py ${h4_all}.structure.e0.txt.lessstructure4:$stable:$down_regulated iCLIP_less_structure:stable:down_regulated ${h4_all}.structure.e0.txt.lessstructure4.KO.stable.down.pdf
# python venn_plot.py ${h4_all}.structure.e0.txt.stablestructure4:$stable:$down_regulated iCLIP_stable_structure:stable:down_regulated ${h4_all}.structure.e0.txt.stablestructure4.KO.stable.down.pdf

# cut -f 1 ${h4_all}.structure.e0.txt.morestructure4 $decay $down_regulated|sort|uniq -c|awk '{if($1==3)print $2}' > ${h4_all}.structure.e0.txt.morestructure4.KO.decay.down.gene.txt
# cut -f 1 ${h4_all}.structure.e0.txt.lessstructure4 $stable $down_regulated|sort|uniq -c|awk '{if($1==3)print $2}' > ${h4_all}.structure.e0.txt.lessstructure4.KO.stable.down.gene.txt
# cut -f 1 ${h4_all}.structure.e0.txt.stablestructure4 $stable $down_regulated|sort|uniq -c|awk '{if($1==3)print $2}' > ${h4_all}.structure.e0.txt.stablestructure4.KO.stable.down.gene.txt


# h6_all=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_all.utr3.bed
# h4_all_minux_h6_tx=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_all_minux_h6_tx.bed
# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' $h6_all $h4_all > $h4_all_minux_h6_tx

# python bed_structure_compare.py $h4_all_minux_h6_tx H4:H6 $shape_sphere:$shape_shield ${h4_all_minux_h6_tx}.structure.e0.txt 0
# python venn_plot.py ${h4_all_minux_h6_tx}.structure.e0.txt.morestructure:$down_regulated:$decay iCLIP_more_structure:down_regulated:decay ${h4_all_minux_h6_tx}.structure.e0.txt.morestructure.KO.decay.png
# python venn_plot.py ${h4_all_minux_h6_tx}.structure.e0.txt.morestructure2:$down_regulated:$decay iCLIP_more_structure:down_regulated:decay ${h4_all_minux_h6_tx}.structure.e0.txt.morestructure2.KO.decay.png
# python venn_plot.py ${h4_all_minux_h6_tx}.structure.e0.txt.lessstructure:$down_regulated:$stable iCLIP_less_structure:down_regulated:stable ${h4_all_minux_h6_tx}.structure.e0.txt.lessstructure.KO.stable.png
# python venn_plot.py ${h4_all_minux_h6_tx}.structure.e0.txt.lessstructure2:$down_regulated:$stable iCLIP_less_structure:down_regulated:stable ${h4_all_minux_h6_tx}.structure.e0.txt.lessstructure2.KO.stable.png
# python venn_plot.py ${h4_all_minux_h6_tx}.structure.e0.txt.stablestructure:$down_regulated:$stable iCLIP_stable_structure:down_regulated:stable ${h4_all_minux_h6_tx}.structure.e0.txt.stablestructure.KO.stable.png
# python venn_plot.py ${h4_all_minux_h6_tx}.structure.e0.txt.stablestructure2:$down_regulated:$stable iCLIP_stable_structure:down_regulated:stable ${h4_all_minux_h6_tx}.structure.e0.txt.stablestructure2.KO.stable.png


# h4_HUR=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-5/bwa/CTK_Procedure/CITS/iCLIP.tag.uniq.clean.CITS.p05.ext20.tranIn.len41.sort.merge.anno.utr3.fimo/fimo.new.HUR_motif2.bed
# h6_HUR=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-6/bwa/CTK_Procedure/CITS/iCLIP.tag.uniq.clean.CITS.p05.ext20.tranIn.len41.sort.merge.anno.utr3.fimo/fimo.new.HUR_motif2.bed
h4_h6_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_HUR_union.bed
# cat $h4_HUR $h6_HUR|awk '{print $1"\t"$2"\t"$3}'|sort|uniq > $h4_h6_HUR
# h4_notiCLIP_HUR=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-5/bwa/CTK_Procedure/CITS/UTR3_minus_iCLIP.fimo/fimo.new.HUR_motif2.bed
h4_h6_HUR_bg=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_HUR_union_bg.bed
# awk 'NR==FNR{a[$1$2$3]=$0;next} !($1$2$3 in a){print $1"\t"$2"\t"$3}' $h4_h6_HUR $h4_notiCLIP_HUR > $h4_h6_HUR_bg
# savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_HUR_union.H4structure.ext20.pdf
# python bed_motif_structure_meta.py $h4_h6_HUR:$h4_h6_HUR_bg H4:bg $shape_sphere:$shape_sphere $savefn 20
# savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_HUR_union.H6structure.ext20.pdf
# python bed_motif_structure_meta.py $h4_h6_HUR:$h4_h6_HUR_bg H6:bg $shape_shield:$shape_shield $savefn 20

# 1st batch, 2 rep
UTR3_bed=/Share2/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.trans.UTR3.bed
tranFa=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.fa
hur_meme=/Share/home/zhangqf7/gongjing/zebrafish/script/zhangting/paris_RBP/motif_CISBP_RNA_narrow/Collapsed.used.HuR2.meme
tranSize=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.size

# h4_12_bed=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/iCLIP.tag.uniq.clean.CITS.p05.ext20.tranIn.len41.sort.merge.anno.utr3.bed
# UTR3_minus_h4_12_bed=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/UTR3_minus_iCLIP.bed
# UTR3_minus_h4_12_fa=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/UTR3_minus_iCLIP.fa
# bedtools subtract -a $UTR3_bed -b $h4_12_bed > $UTR3_minus_h4_12_bed
# bedtools getfasta -fi $tranFa -bed $UTR3_minus_h4_12_bed -fo $UTR3_minus_h4_12_fa

# shuffle_bed=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/UTR3_minus_iCLIP.bg/UTR3_minus_iCLIP.bed
# shuffle_fa=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/UTR3_minus_iCLIP.bg/UTR3_minus_iCLIP.fa
# shuffle_bg=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/UTR3_minus_iCLIP.bg/UTR3_minus_iCLIP.bg
# shuffleBed -seed 1234 -i $UTR3_minus_h4_12_bed -chrom -g $tranSize > $shuffle_bed
# bedtools getfasta -fi $tranFa -bed $shuffle_bed -fo $shuffle_fa
# fasta-get-markov $shuffle_fa $shuffle_bg
# fimo_dir=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/UTR3_minus_iCLIP.fimo
# fimo --oc $fimo_dir --thresh 0.001 --norc --bgfile $shuffle_bg $hur_meme $UTR3_minus_h4_12_fa


# h4_12_bed=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h6_full_12_all/CITS/iCLIP.tag.uniq.clean.CITS.p05.ext20.tranIn.len41.sort.merge.anno.utr3.bed
# UTR3_minus_h4_12_bed=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h6_full_12_all/CITS/UTR3_minus_iCLIP.bed
# UTR3_minus_h4_12_fa=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h6_full_12_all/CITS/UTR3_minus_iCLIP.fa
# bedtools subtract -a $UTR3_bed -b $h4_12_bed > $UTR3_minus_h4_12_bed
# bedtools getfasta -fi $tranFa -bed $UTR3_minus_h4_12_bed -fo $UTR3_minus_h4_12_fa

# shuffle_bed=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h6_full_12_all/CITS/UTR3_minus_iCLIP.bg/UTR3_minus_iCLIP.bed
# shuffle_fa=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h6_full_12_all/CITS/UTR3_minus_iCLIP.bg/UTR3_minus_iCLIP.fa
# shuffle_bg=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h6_full_12_all/CITS/UTR3_minus_iCLIP.bg/UTR3_minus_iCLIP.bg
# shuffleBed -seed 1234 -i $UTR3_minus_h4_12_bed -chrom -g $tranSize > $shuffle_bed
# bedtools getfasta -fi $tranFa -bed $shuffle_bed -fo $shuffle_fa
# fasta-get-markov $shuffle_fa $shuffle_bg
# fimo_dir=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h6_full_12_all/CITS/UTR3_minus_iCLIP.fimo
# fimo --oc $fimo_dir --thresh 0.001 --norc --bgfile $shuffle_bg $hur_meme $UTR3_minus_h4_12_fa


# h4_12_HUR=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/iCLIP.tag.uniq.clean.CITS.p05.c6.ext20.tranIn.len41.sort.merge.anno.utr3.fimo2/fimo.new.txt
# h6_12_HUR=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h6_full_12_all/CITS/iCLIP.tag.uniq.clean.CITS.p05.c6.ext20.tranIn.len41.sort.merge.anno.utr3.fimo2/fimo.new.txt
# h4_h6_12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_12_HUR_union.bed
# cat $h4_12_HUR $h6_12_HUR|awk '{print $1"\t"$2"\t"$3}'|sort|uniq > $h4_h6_12_HUR
# h4_12_notiCLIP_HUR=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/iCLIP.tag.uniq.clean.CITS.p05.c8.ext20.tranIn.len41.sort.merge.anno.utr3.UTR3minusiCLIP.fimo2/fimo.new.txt
# h4_h6_12_HUR_bg=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_12_HUR_union_bg.bed
# awk 'NR==FNR{a[$1$2$3]=$0;next} !($1$2$3 in a){print $1"\t"$2"\t"$3}' $h4_h6_12_HUR $h4_12_notiCLIP_HUR > $h4_h6_12_HUR_bg
# savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_12_HUR_union.H4structure.ext20.c6.pdf

# python bed_motif_structure_meta.py $h4_h6_12_HUR:$h4_h6_12_HUR_bg H4:bg $shape_sphere:$shape_sphere $savefn 20


# h4_h6_12_HUR_shuffle=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_12_HUR_union.shuffle.bed
# shuffleBed -seed 1234 -i $h4_h6_12_HUR -chrom -g $tranSize > $h4_h6_12_HUR_shuffle

# h4_h6_HUR_in_h4_h6_12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_12_HUR_union.inSingle.bed
# awk 'NR==FNR{a[$1$2$3]=$0;next} $1$2$3 in a{print $1"\t"$2"\t"$3}' $h4_h6_HUR $h4_h6_12_HUR > $h4_h6_HUR_in_h4_h6_12_HUR
# h4_h6_HUR_notin_h4_h6_12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_12_HUR_union.notinSingle.bed
# awk 'NR==FNR{a[$1$2$3]=$0;next} !($1$2$3 in a){print $1"\t"$2"\t"$3}' $h4_h6_HUR $h4_h6_12_HUR > $h4_h6_HUR_notin_h4_h6_12_HUR

# h4_h6_HUR_bg_in_new
# awk 'NR==FNR{a[$1$2$3]=$0;next} $1$2$3 in a{print $1"\t"$2"\t"$3}' $h4_h6_12_HUR_bg $h4_h6_HUR_bg > 

# python bed_motif_structure_meta.py $h4_h6_12_HUR:$h4_h6_12_HUR_bg:$h4_h6_12_HUR_shuffle:$h4_h6_HUR_in_h4_h6_12_HUR:$h4_h6_HUR_notin_h4_h6_12_HUR H4:bg:shuffle:H4_in_before:H4_notin_before $shape_sphere:$shape_sphere:$shape_sphere:$shape_sphere:$shape_sphere $savefn 20


h4_all_rep12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_all_rep12.c6.utr3.HuRmotif.bed
h6_all_rep12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_all_rep12.c6.utr3.HuRmotif.bed
h4_h6_rep12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_h6_all_rep12.c6.utr3.HuRmotif.bed
# cat $h4_all_rep12_HUR $h6_all_rep12_HUR|awk '{print $1"\t"$2"\t"$3}'|sort|uniq > $h4_h6_rep12_HUR
# awk 'NR==FNR{a[$1]=$0;next} !($1 in a){print $1"\t"$2"\t"$3}' $decay $h4_h6_rep12_HUR > ${h4_h6_rep12_HUR}.decay
# awk 'NR==FNR{a[$1]=$0;next} !($1 in a){print $1"\t"$2"\t"$3}' $stable $h4_h6_rep12_HUR > ${h4_h6_rep12_HUR}.stable
savefn=${h4_h6_rep12_HUR}.decay_stable_diff.pdf
# python bed_motif_structure_meta.py ${h4_h6_rep12_HUR}.decay:${h4_h6_rep12_HUR}.decay:${h4_h6_rep12_HUR}.stable:${h4_h6_rep12_HUR}.stable decay:decay:stable:stable $shape_shield:$shape_sphere:$shape_shield:$shape_sphere $savefn 20

# awk 'NR==FNR{a[$1]=$0;next} !($1 in a){print $1"\t"$2"\t"$3}' ${h4_h6_rep12_HUR}.decay ${h4_h6_rep12_HUR}.stable > ${h4_h6_rep12_HUR}.stable.nocommon
# awk 'NR==FNR{a[$1]=$0;next} !($1 in a){print $1"\t"$2"\t"$3}' ${h4_h6_rep12_HUR}.stable ${h4_h6_rep12_HUR}.decay > ${h4_h6_rep12_HUR}.decay.nocommon

savefn=${h4_h6_rep12_HUR}.decay_stable_diff.nocomon.pdf
# python bed_motif_structure_meta.py ${h4_h6_rep12_HUR}.decay.nocommon:${h4_h6_rep12_HUR}.decay.nocommon:${h4_h6_rep12_HUR}.stable.nocommon:${h4_h6_rep12_HUR}.stable.nocommon decay:decay:stable:stable $shape_shield:$shape_sphere:$shape_shield:$shape_sphere $savefn 20

# h4_minus_h6_rep12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_minus_h6_all_rep12.c6.utr3.HuRmotif.bed
# awk 'NR==FNR{a[$1$2$3]=$0;next} !($1$2$3 in a){print $1"\t"$2"\t"$3}' $h6_all_rep12_HUR $h4_all_rep12_HUR > $h4_minus_h6_rep12_HUR
# awk 'NR==FNR{a[$1]=$0;next} !($1 in a){print $1"\t"$2"\t"$3}' $decay $h4_minus_h6_rep12_HUR > ${h4_minus_h6_rep12_HUR}.decay
# awk 'NR==FNR{a[$1]=$0;next} !($1 in a){print $1"\t"$2"\t"$3}' $stable $h4_minus_h6_rep12_HUR > ${h4_minus_h6_rep12_HUR}.stable
# savefn=${h4_minus_h6_rep12_HUR}.decay_stable_diff.pdf
# python bed_motif_structure_meta.py ${h4_minus_h6_rep12_HUR}.decay:${h4_minus_h6_rep12_HUR}.decay:${h4_minus_h6_rep12_HUR}.stable:${h4_minus_h6_rep12_HUR}.stable decay:decay:stable:stable $shape_shield:$shape_sphere:$shape_shield:$shape_sphere $savefn 20


m6A=/Share2/home/zhangqf7/gongjing/zebrafish/result/decay_sets/DF2_targets_up1.5.maternal.ID
mir430=/Share2/home/zhangqf7/gongjing/zebrafish/result/decay_sets/mir430.maternal.targetsID
elavl1a=/Share2/home/zhangqf7/gongjing/zebrafish/result/decay_sets/h4_all_rep12.c6.utr3.bed.structure.e0.txt.morestructure4.decay.lessstructure4.stable.ensembl.txt
non_optimal=/Share2/home/zhangqf7/gongjing/zebrafish/result/decay_sets/cai_res_bottom10percent_Ensemblgene.txt
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/decay_sets/overlap.pdf
# python venn_plot.py $m6A:$mir430:$elavl1a:$non_optimal m6A:mir430:elavl1a:non_optimal $savefn

amanitin=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq0117/Amanitin_H2_H6.DE.decay.lfc1.2.q0.05.genes.txt
morestructure4=/Share/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_all_rep12.c6.utr3.bed.structure.e0.txt.morestructure4
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq0117/Amanitin_H2_H6.DE.decay.lfc1.2.q0.05.genes.iclipmorestructure4.overlap.pdf
python venn_plot.py $amanitin:$morestructure4:$decay Maternal_dependent_decay:More_structure:Maternal_decay $savefn


