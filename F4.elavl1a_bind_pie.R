fimo=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/abs/sphere_shield/fimo_out-utr3/fimo.txt
grep 'HUR_motif2' $fimo|cut -f3|sed 's/:/\t/g' |cut -f1|sort -u|wc #813

cut -f1 /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/abs/sphere_shield/window-anno_utr3.bed|sort -u|wc #2635

fimo=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_005_new/abs/sphere_shield/fimo_out-utr3/fimo.txt
grep 'HUR_motif2' $fimo|cut -f3|sed 's/:/\t/g' |cut -f1|sort -u|wc #1423
cut -f1 /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_005_new/abs/sphere_shield/window-anno_utr3.bed|sort -u|wc #3667


fimo=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_001_new/abs/sphere_shield/fimo_out-utr3/fimo.txt
grep 'HUR_motif2' $fimo|cut -f3|sed 's/:/\t/g' |cut -f1|sort -u|wc #1024
cut -f1 /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_001_new/abs/sphere_shield/window-anno_utr3.bed|sort -u|wc #2960


fimo=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_005_new/abs/sphere_shield/fimo_out-utr3/fimo.txt
grep 'HUR_motif2' $fimo|cut -f3|sed 's/:/\t/g' |cut -f1|sort -u|wc #1895
cut -f1 /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_005_new/abs/sphere_shield/window-anno_utr3.bed|sort -u|wc #4023


[zhangqf7@ZIO01 sphere_shield]$ pwd
/Share2/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/sphere_shield
[zhangqf7@ZIO01 window-anno_utr3]$ grep 'HUR_motif2' fimo.txt|cut -f3|sed 's/:/\t/g' |cut -f1|sort -u|wc
   1895    1895   21845
[zhangqf7@ZIO01 sphere_shield]$ cut -f1 window-anno_utr3.bed|sort -u|wc
   4023    4023   47079

setwd('/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/pic/FigureS4E')
png('FigureS4E-01_001.png')
par(lwd=3)
pie(c(813,2635), col=c('red','lightgrey'),labels=NA)
dev.off()


png('FigureS4E-01_005.png')
par(lwd=3)
pie(c(1423,3667), col=c('red','lightgrey'),labels=NA)
dev.off()


png('FigureS4E-005_001.png')
par(lwd=3)
pie(c(1024,2960), col=c('red','lightgrey'),labels=NA)
dev.off()


png('FigureS4E-005_005.png')
par(lwd=3)
pie(c(1895,4023), col=c('red','lightgrey'),labels=NA)
dev.off()