h4_12_HUR=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/iCLIP.tag.uniq.clean.CITS.p05.c6.ext20.tranIn.len41.sort.merge.anno.utr3.fimo2/fimo.new.txt
h6_12_HUR=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h6_full_12_all/CITS/iCLIP.tag.uniq.clean.CITS.p05.c6.ext20.tranIn.len41.sort.merge.anno.utr3.fimo2/fimo.new.txt
h4_h6_12_HUR=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_12_HUR_union.bed
cat $h4_12_HUR $h6_12_HUR|awk '{print $1"\t"$2"\t"$3}'|sort|uniq > $h4_h6_12_HUR
h4_12_notiCLIP_HUR=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/h4_full_12_all/CITS/iCLIP.tag.uniq.clean.CITS.p05.c8.ext20.tranIn.len41.sort.merge.anno.utr3.UTR3minusiCLIP.fimo2/fimo.new.txt
h4_h6_12_HUR_bg=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_12_HUR_union_bg.bed
awk 'NR==FNR{a[$1$2$3]=$0;next} !($1$2$3 in a){print $1"\t"$2"\t"$3}' $h4_h6_12_HUR $h4_12_notiCLIP_HUR > $h4_h6_12_HUR_bg
savefn=/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/H4_H6_12_HUR_union.H4structure.ext20.c6.pdf

python bed_motif_structure_meta.py $h4_h6_12_HUR:$h4_h6_12_HUR_bg H4:bg $shape_sphere:$shape_sphere $savefn 20