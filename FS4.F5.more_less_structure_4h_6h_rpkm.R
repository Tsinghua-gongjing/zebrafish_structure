 # setwd('/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/pic/Figure2G-H')
# path="/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/pic/Figure2D/01_005-"
# output="01_005-"
# abs.sphere <- read.table(paste(path,'res_abs/HUR_motif2-sphere.out',sep=''))
# abs.shield <- read.table(paste(path,'res_abs/HUR_motif2-sshield.out',sep=''))

# x1 <- 0:40
# png(paste(output,"hur-abs-minus.png",sep=""),width=1000,height=500)
# par(cex.axis=2,mgp=c(3,3,0),lwd=4,cex.lab=2,mar=c(8,9,2,2))
# plot(-1,type="l",lwd=3,col="black",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,40),ylim=c(-0.1,0.1));
# lines(x1,abs.shield-abs.sphere,lwd=5,col="black",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
# points(x1,abs.shield-abs.sphere,col="black",cex=1)
# abline(v=17,lwd=2, lty=2)
# abline(v=23,lwd=2, lty=2)
# axis(1,at=c(0,20,40),labels=c(-20,"Site",20),tcl=-2, las=1,lwd=4, pos=-0.1)
# axis(2,at=c(-0.1,-0.05,0,0.05,0.1),labels=c(-0.1,-0.05,0,0.05,0.1),tcl=-2, las=1,lwd=4, pos=0)
# dev.off()


# up.sphere <- read.table(paste(path,'res_up/HUR_motif2-sphere.out',sep=''))
# up.shield <- read.table(paste(path,'res_up/HUR_motif2-sshield.out',sep=''))

# png(paste(output,"hur-up-minus.png",sep=""),width=1000,height=500)
# par(cex.axis=2,mgp=c(3,3,0),lwd=4,cex.lab=2,mar=c(8,9,2,2))
# plot(-1,type="l",lwd=3,col="firebrick",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,40),ylim=c(-0.1,0.1));
# lines(x1,up.shield-up.sphere,lwd=5,col="firebrick",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
# points(x1,up.shield-up.sphere,col="firebrick",cex=1)
# abline(v=17,lwd=2, lty=2)
# abline(v=23,lwd=2, lty=2)
# axis(1,at=c(0,20,40),labels=c(-20,"Site",20),tcl=-2, las=1,lwd=4, pos=-0.1)
# axis(2,at=c(-0.1,-0.05,0,0.05,0.1),labels=c(-0.1,-0.05,0.00,0.05,0.1),tcl=-2, las=1,lwd=4, pos=0)
# dev.off()


# down.sphere <- read.table(paste(path,'res_down/HUR_motif2-sphere.out',sep=''))
# down.shield <- read.table(paste(path,'res_down/HUR_motif2-sshield.out',sep=''))

# png(paste(output,"hur-down-minus.png",sep=""),width=1000,height=500)
# par(cex.axis=2,mgp=c(3,3,0),lwd=4,cex.lab=2,mar=c(8,9,2,2))
# plot(-1,type="l",lwd=5,col="darkblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,40),ylim=c(-0.1,0.1));
# lines(x1,down.shield-down.sphere,lwd=5,col="darkblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
# points(x1,down.shield-down.sphere,col="darkblue",cex=1)
# abline(v=17,lwd=2, lty=2)
# abline(v=23,lwd=2, lty=2)
# axis(1,at=c(0,20,40),labels=c(-20,"Site",20),tcl=-2, las=1,lwd=4, pos=-0.1)
# axis(2,at=c(-0.1,-0.05,0,0.05,0.1),labels=c(-0.1,-0.05,0.00,0.05,0.1),tcl=-2, las=1,lwd=4, pos=0)
# dev.off()





#######################################基因集合的比较
setwd('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/Figure-structure-main/Figure2G-H')

# use iclip (old)
# hur.trans = read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/iCLIP/overlap/shi-zi-5.shi-zi-6.ol.structure.bed',stringsAsFactors=F,sep='\t',header=F)
# hur.up.trans <- hur.trans[which(hur.trans$V8=="up"), "V1"]
# hur.down.trans <- hur.trans[which(hur.trans$V8=="down"), "V1"]
# print(length(hur.up.trans))
# print(length(hur.down.trans))

# use iclip, new rep1/2
hur.up.trans <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/h4/h4_all_rep12.c6.utr3.bed.structure.e0.txt.lessstructure4',stringsAsFactors=F,sep='\t',header=F)
hur.down.trans <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/h4/h4_all_rep12.c6.utr3.bed.structure.e0.txt.morestructure4',stringsAsFactors=F,sep='\t',header=F)
hur.stable.trans <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/h4/h4_all_rep12.c6.utr3.bed.structure.e0.txt.stablestructure4',stringsAsFactors=F,sep='\t',header=F)
hur.up.trans <- hur.up.trans$V1
hur.down.trans <- hur.down.trans$V1
hur.stable.trans <- hur.stable.trans$V1

# use predict
# hur.up.trans <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/enrich_new/new.sphere_shield.up-enrich.utr3.txt',stringsAsFactors=F,sep='\t',row.names=1,header=T)
# hur.down.trans <- read.table('/Share2/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/enrich_new/new.sphere_shield.down-enrich.utr3.txt',stringsAsFactors=F,sep='\t',row.names=1,header=T)
# # hur.up.trans <- unlist(strsplit(hur.up.trans['HUR_motif2',4],"; "))
# # hur.down.trans <- unlist(strsplit(hur.down.trans['HUR_motif2',4],"; "))
# hur.up.trans <- unlist(strsplit(strsplit(hur.up.trans['HUR_motif2','transNum'], "; ")[[1]][2], ", "))
# hur.down.trans <- unlist(strsplit(strsplit(hur.down.trans['HUR_motif2','transNum'], "; ")[[1]][2], ", "))
# print(length(hur.up.trans))
# print(length(hur.down.trans))


hur.up.notdown.trans = setdiff(hur.up.trans, hur.down.trans)
hur.down.notup.trans = setdiff(hur.down.trans, hur.up.trans)
print(length(hur.up.notdown.trans))
print(length(hur.down.notup.trans))
print(length(hur.stable.trans))

output <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq_ctrl_h46/output_score.txt',sep='\t',stringsAsFactors=F,header=T,row.names=1)
# output <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-hur-mo-h46/output_score.txt',sep='\t',stringsAsFactors=F,header=T,row.names=1)
output <- output[complete.cases(output),]
fc <- output[,'log2.Fold_change..normalized']
names(fc) <- rownames(output)



# pdf(paste("hur_utr3_01_005-cv.pdf",sep=""),width=10)
pdf(paste("hur_utr3_005_005-cv.iclip.nocommon.pdf",sep=""), width=10)
# pdf(paste("hur_utr3_01_001-cv.pdf",sep=""), width=10)
par(cex.lab=1.5, cex.axis=1.5, mgp=c(3.5,1,0), mar=c(5,6,5,5))
# up <- fc[hur.up.trans]
up <- fc[hur.up.notdown.trans]
# down <- fc[hur.down.trans]
down <- fc[hur.down.notup.trans]
# change <- fc[unique(c(hur.up.trans,hur.down.trans))]
change <- fc[hur.stable.trans] # stable
stable <- fc[hur.stable.trans] # stable

x1 <- ecdf(up + seq(0.00000001,by=0.000000001,length.out=length(up)))
x2 <- ecdf(down+ seq(0.00000001, by=0.000000001,length.out=length(down)))
x3 <- ecdf(change + seq(0.00000001, by=0.000000001,length.out=length(change)))

plot(x2,lwd=1,col="white",main="",xlab="log2(6h/4h)",ylab="Culmulative Frequency", yaxt="n",xlim=c(-2,2), xaxt="n", pch=19, cex=1.5)
axis(1,at=seq(-2,2,1),labels=rep("",5),tcl=-1.2,lwd=3)
axis(2,at=seq(0,1,0.2),labels=rep("",6),tcl=-1.2,lwd=3)
lines(x1, col="green",do.points=F, lwd=3)
lines(x2, col="firebrick",do.points=F, lwd=3)
lines(x3, col="darkblue",do.points=F, lwd=3)
box(lwd=3)
dev.off()

print("up vs down")
wilcox.test(up, down)$p.value
wilcox.test(up, down, alternative='less')$p.value
wilcox.test(up, down, alternative='greater')$p.value
# wilcox.test(up, stable)$p.value
# wilcox.test(down, change )$p.value
# ks.test(down, stable )$p.value
print("up vs stable")
wilcox.test(up, stable)$p.value
wilcox.test(up, stable, alternative='less')$p.value
wilcox.test(up, stable, alternative='greater')$p.value

print("down vs stable")
wilcox.test(down, stable)$p.value
wilcox.test(down, stable, alternative='less')$p.value
wilcox.test(down, stable, alternative='greater')$p.value


pdf(paste("hur_utr3_005_005-cv.iclip.nocommon.box.pdf",sep=""))
boxplot(up, down, change, names=c('less structure', 'more structure', 'stable'), col = c("green","firebrick","darkblue"))
dev.off()









































