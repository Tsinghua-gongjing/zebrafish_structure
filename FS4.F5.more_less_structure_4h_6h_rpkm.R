#######################################基因集合的比较
setwd('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/Figure-structure-main/Figure2G-H')

# use iclip (old)
# hur.trans = read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/iCLIP/overlap/shi-zi-5.shi-zi-6.ol.structure.bed',stringsAsFactors=F,sep='\t',header=F)
# hur.up.trans <- hur.trans[which(hur.trans$V8=="up"), "V1"]
# hur.down.trans <- hur.trans[which(hur.trans$V8=="down"), "V1"]
# print(length(hur.up.trans))
# print(length(hur.down.trans))

# use iclip, new rep1/2
# hur.up.trans <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/h4/h4_all_rep12.c6.utr3.bed.structure.e0.txt.lessstructure4',stringsAsFactors=F,sep='\t',header=F)
# hur.down.trans <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/h4/h4_all_rep12.c6.utr3.bed.structure.e0.txt.morestructure4',stringsAsFactors=F,sep='\t',header=F)
# # hur.stable.trans <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/h4/h4_all_rep12.c6.utr3.bed.structure.e0.txt.stablestructure4',stringsAsFactors=F,sep='\t',header=F)
# hur.up.trans <- hur.up.trans$V1
# hur.down.trans <- hur.down.trans$V1
# # hur.stable.trans <- hur.stable.trans$V1

# use predict
hur.up.trans <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/enrich_new/new.sphere_shield.up-enrich.utr3.txt',stringsAsFactors=F,sep='\t',row.names=1,header=T)
hur.down.trans <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/enrich_new/new.sphere_shield.down-enrich.utr3.txt',stringsAsFactors=F,sep='\t',row.names=1,header=T)
# hur.up.trans <- unlist(strsplit(hur.up.trans['HUR_motif2',4],"; "))
# hur.down.trans <- unlist(strsplit(hur.down.trans['HUR_motif2',4],"; "))
hur.up.trans <- unlist(strsplit(strsplit(hur.up.trans['HUR_motif2','transNum'], "; ")[[1]][2], ", "))
hur.down.trans <- unlist(strsplit(strsplit(hur.down.trans['HUR_motif2','transNum'], "; ")[[1]][2], ", "))
print(length(hur.up.trans))
print(length(hur.down.trans))


hur.up.notdown.trans = setdiff(hur.up.trans, hur.down.trans)
hur.down.notup.trans = setdiff(hur.down.trans, hur.up.trans)
print(length(hur.up.notdown.trans))
print(length(hur.down.notup.trans))

# output <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq_ctrl_h46/output_score.txt',sep='\t',stringsAsFactors=F,header=T,row.names=1)
output <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/RNAseq/DEGresult/control_4h_vs_6h/output_score.txt',sep='\t',stringsAsFactors=F,header=T,row.names=1)
output <- output[complete.cases(output),]
fc <- output[,'log2.Fold_change..normalized']
names(fc) <- rownames(output)


# pdf(paste("/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/F5.hur_utr3_005_005-cv.iclip.nocommon.pdf",sep=""), width=10)
pdf(paste("/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/F5.hur_utr3_005_005-cv.usepredict.nocommon.pdf",sep=""), width=10)
par(cex.lab=1.5, cex.axis=1.5, mgp=c(3.5,1,0), mar=c(5,6,5,5))
# up <- fc[hur.up.trans]
up <- fc[hur.up.notdown.trans]
# down <- fc[hur.down.trans]
down <- fc[hur.down.notup.trans]
# change <- fc[unique(c(hur.up.trans,hur.down.trans))]
# change <- fc[hur.stable.trans] # stable
# stable <- fc[hur.stable.trans] # stable

x1 <- ecdf(up + seq(0.00000001,by=0.000000001,length.out=length(up)))
x2 <- ecdf(down+ seq(0.00000001, by=0.000000001,length.out=length(down)))
# x3 <- ecdf(change + seq(0.00000001, by=0.000000001,length.out=length(change)))

plot(x2,lwd=1,col="white",main="",xlab="log2(6h/4h)",ylab="Culmulative Frequency", yaxt="n",xlim=c(-2,2), xaxt="n", pch=19, cex=1.5)
axis(1,at=seq(-2,2,1),labels=rep("",5),tcl=-1.2,lwd=3)
axis(2,at=seq(0,1,0.2),labels=rep("",6),tcl=-1.2,lwd=3)
lines(x1, col="green",do.points=F, lwd=3)
lines(x2, col="firebrick",do.points=F, lwd=3)
# lines(x3, col="darkblue",do.points=F, lwd=3)
box(lwd=3)
dev.off()

print("up vs down")
wilcox.test(up, down)$p.value
wilcox.test(up, down, alternative='less')$p.value
wilcox.test(up, down, alternative='greater')$p.value

# print("up vs stable")
# wilcox.test(up, stable)$p.value
# wilcox.test(up, stable, alternative='less')$p.value
# wilcox.test(up, stable, alternative='greater')$p.value

# print("down vs stable")
# wilcox.test(down, stable)$p.value
# wilcox.test(down, stable, alternative='less')$p.value
# wilcox.test(down, stable, alternative='greater')$p.value


# iclip
# [zhangqf7@ZIO01 zebrafish_structure]$ Rscript FS4.F5.more_less_structure_4h_6h_rpkm.R
# [1] 248
# [1] 817
# [1] 1613
# null device
#           1
# [1] "up vs down"
# [1] 6.897413e-15
# [1] 1
# [1] 3.448707e-15
# [1] "up vs stable"
# [1] 1.798221e-07
# [1] 0.9999999
# [1] 8.991103e-08
# [1] "down vs stable"
# [1] 6.644164e-08
# [1] 3.322082e-08
# [1] 1


# # use predict
# [zhangqf7@ZIO01 zebrafish_structure]$ Rscript FS4.F5.more_less_structure_4h_6h_rpkm.R
# [1] 870
# [1] 1485
# [1] 409
# [1] 1024
# null device
#           1
# [1] "up vs down"
# [1] 3.948819e-17
# [1] 1
# [1] 1.97441e-17

