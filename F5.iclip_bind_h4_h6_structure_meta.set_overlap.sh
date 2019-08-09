h4_all=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/h4/h4_all_rep12.c6.utr3.bed
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

decay=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_New_refseq/decay_id.txt
stable=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_New_refseq/stable_id.txt
zygotic=/Share2/home/zhangqf7/gongjing/zebrafish/result/RNAseq_New_refseq/zygotic_id.txt

# python venn_plot.py ${h4_all}.structure.e0.txt.morestructure4:$decay:$stable iCLIP_more_structure:decay:stable ${h4_all}.structure.e0.txt.morestructure4.KO.decay.pdf
# python venn_plot.py ${h4_all}.structure.e0.txt.lessstructure4:$decay:$stable iCLIP_less_structure:decay:stable ${h4_all}.structure.e0.txt.lessstructure4.KO.stable.pdf
# python venn_plot.py ${h4_all}.structure.e0.txt.stablestructure4:$decay:$stable iCLIP_stable_structure:decay:stable ${h4_all}.structure.e0.txt.stablestructure4.KO.stable.pdf

# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $1}' $decay ${h4_all}.structure.e0.txt.morestructure4 > ${h4_all}.structure.e0.txt.morestructure4.decay
# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' ${h4_all}.structure.e0.txt.morestructure4.decay ${h4_all}|awk '{if($3-$2==41)print}' > ${h4_all}.structure.e0.txt.morestructure4.decay.iclip.bed 

# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $1}' $stable ${h4_all}.structure.e0.txt.lessstructure4 > ${h4_all}.structure.e0.txt.lessstructure4.stable
# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' ${h4_all}.structure.e0.txt.lessstructure4.stable ${h4_all}|awk '{if($3-$2==41)print}' > ${h4_all}.structure.e0.txt.lessstructure4.stable.iclip.bed 

# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $1}' $stable ${h4_all}.structure.e0.txt.stablestructure4 > ${h4_all}.structure.e0.txt.stablestructure4.stable
# awk 'NR==FNR{a[$1]=$0;next} $1 in a{print $0}' ${h4_all}.structure.e0.txt.stablestructure4.stable ${h4_all}|awk '{if($3-$2==41)print}' > ${h4_all}.structure.e0.txt.stablestructure4.stable.iclip.bed 

ol1=${h4_all}.structure.e0.txt.morestructure4.decay.iclip.bed
ol2=${h4_all}.structure.e0.txt.lessstructure4.stable.iclip.bed
ol3=${h4_all}.structure.e0.txt.stablestructure4.stable.iclip.bed 


shape_sphere=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/icSHAPE/shape/sphere.icshape.w200.s30.T2.t200.out
shape_shield=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/icSHAPE/shape/shield.icshape.w200.s30.T2.t200.out

python bed_motif_structure_meta.py $ol1:$ol1 H4:H6 $shape_sphere:$shape_shield ${ol1}.structure.ext0.pdf 0 
python bed_motif_structure_meta.py $ol2:$ol2 H4:H6 $shape_sphere:$shape_shield ${ol2}.structure.ext0.pdf 0 
python bed_motif_structure_meta.py $ol3:$ol3 H4:H6 $shape_sphere:$shape_shield ${ol3}.structure.ext0.pdf 0 