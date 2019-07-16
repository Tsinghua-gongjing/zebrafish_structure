shape_sphere=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/icSHAPE/shape/sphere.icshape.w200.s30.T2.t200.out
shape_shield=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/icSHAPE/shape/shield.icshape.w200.s30.T2.t200.out

h4_12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/h4/fimo/fimo.new.txt
h6_12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/h6/fimo/fimo.new.txt
h4_h6_12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/H4_H6_12_HUR_union.bed
cat $h4_12_HUR $h6_12_HUR|awk '{print $1"\t"$2"\t"$3}'|sort|uniq > $h4_h6_12_HUR
h4_12_notiCLIP_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/h4/notiCLIP_fimo/fimo.new.txt
h4_h6_12_HUR_bg=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/H4_H6_12_HUR_union_bg.bed
awk 'NR==FNR{a[$1$2$3]=$0;next} !($1$2$3 in a){print $1"\t"$2"\t"$3}' $h4_h6_12_HUR $h4_12_notiCLIP_HUR > $h4_h6_12_HUR_bg
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/H4_H6_12_HUR_union.H4structure.ext20.c6.pdf

python bed_motif_structure_meta.py $h4_h6_12_HUR:$h4_h6_12_HUR_bg H4:bg $shape_sphere:$shape_sphere $savefn 20 # bed_motif_structure_meta: bed_meta2