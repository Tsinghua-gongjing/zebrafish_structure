setwd('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/Figure-structure-main/Figure3C')


ctrl_h2.rep1 <- read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/control_h2-rep1-rpkm.txt",sep='\t',row.names=1,header=F)
ctrl_h2.rep2 <- read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/control_h2-rep2-rpkm.txt",sep='\t',row.names=1,header=F)

ctrl_h4.rep1 <- read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/control_h4-rep1-rpkm.txt",sep='\t',stringsAsFactors=F,header=F,row.names=1)
ctrl_h4.rep2 <- read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/control_h4-rep2-rpkm.txt",sep='\t',stringsAsFactors=F,header=F,row.names=1)

ctrl_h6.rep1 <- read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/control_h6-rep1-rpkm.txt",sep='\t',row.names=1,header=F)
ctrl_h6.rep2 <- read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/control_h6-rep2-rpkm.txt",sep='\t',row.names=1,header=F)

hur_mo_h2.rep1 <- read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/hur_MO-h2_rep1-rpkm.txt",sep='\t',row.names=1,header=F)
hur_mo_h2.rep2 <- read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/hur_MO-h2_rep2-rpkm.txt",sep='\t',row.names=1,stringsAsFactors=F)

hur_mo_h4.rep1 <-  read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/hur_MO-h4_rep1-rpkm.txt",sep='\t',stringsAsFactors=F,header=F,row.names=1)
hur_mo_h4.rep2 <-  read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/hur_MO-h4_rep2-rpkm.txt",sep='\t',stringsAsFactors=F,header=F,row.names=1)

hur_mo_h6.rep1 <- read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/hur-MO-h6_rep1f-rpkm.txt",sep='\t',row.names=1,header=F)
hur_mo_h6.rep2 <- read.table("/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-mapping-quantifing-DEGs文件/hur-MO-h6_rep2f-rpkm.txt",sep='\t',row.names=1,header=F)

name <- unique(c(rownames(ctrl_h2.rep1),rownames(ctrl_h2.rep2),rownames(ctrl_h4.rep1),rownames(ctrl_h4.rep2),rownames(ctrl_h6.rep1),rownames(ctrl_h6.rep2),rownames(hur_mo_h2.rep1),rownames(hur_mo_h2.rep2),rownames(hur_mo_h4.rep1),rownames(hur_mo_h4.rep2),rownames(hur_mo_h6.rep1),rownames(hur_mo_h6.rep2)))


#################count
count <- matrix(0, nrow=length(name), ncol=12)
rownames(count) <- name
count[rownames(ctrl_h2.rep1),1] <- ctrl_h2.rep1[,2]+ctrl_h2.rep1[,3]
count[rownames(ctrl_h2.rep2),2] <- ctrl_h2.rep2[,2]+ctrl_h2.rep2[,3]
count[rownames(ctrl_h4.rep1),3] <- ctrl_h4.rep1[,2]+ctrl_h4.rep1[,3]
count[rownames(ctrl_h4.rep2),4] <- ctrl_h4.rep2[,2]+ctrl_h4.rep2[,3]
count[rownames(ctrl_h6.rep1),5] <- ctrl_h6.rep1[,2]+ctrl_h6.rep1[,3]
count[rownames(ctrl_h6.rep2),6] <- ctrl_h6.rep2[,2]+ctrl_h6.rep2[,3]

count[rownames(hur_mo_h2.rep1),7] <- hur_mo_h2.rep1[,2]+hur_mo_h2.rep1[,3]
count[rownames(hur_mo_h2.rep2),8] <- hur_mo_h2.rep2[,2]+hur_mo_h2.rep2[,3]
count[rownames(hur_mo_h4.rep1),9] <- hur_mo_h4.rep1[,2]+hur_mo_h4.rep1[,3]
count[rownames(hur_mo_h4.rep2),10] <- hur_mo_h4.rep2[,2]+hur_mo_h4.rep2[,3]
count[rownames(hur_mo_h6.rep1),11] <- hur_mo_h6.rep1[,2]+hur_mo_h6.rep1[,3]
count[rownames(hur_mo_h6.rep2),12] <- hur_mo_h6.rep2[,2]+hur_mo_h6.rep2[,3]
colnames(count) <- c('ctrl-h2_rep1','ctrl-h2_rep2','ctrl-h4_rep1','ctrl-h4_rep2','ctrl-h6_rep1','ctrl-h6_rep2','hur-h2_rep1','hur-h2_rep2','hur-h4_rep1','hur-h4_rep2','hur-h6_rep1','hur-h6_rep2')

#################rpkm
rpkm <- matrix(0, nrow=length(name), ncol=12)
rownames(rpkm) <- name
rpkm[rownames(ctrl_h2.rep1),1] <- ctrl_h2.rep1[,4]
rpkm[rownames(ctrl_h2.rep2),2] <- ctrl_h2.rep2[,4]
rpkm[rownames(ctrl_h4.rep1),3] <- ctrl_h4.rep1[,4]
rpkm[rownames(ctrl_h4.rep2),4] <- ctrl_h4.rep2[,4]
rpkm[rownames(ctrl_h6.rep1),5] <- ctrl_h6.rep1[,4]
rpkm[rownames(ctrl_h6.rep2),6] <- ctrl_h6.rep2[,4]

rpkm[rownames(hur_mo_h2.rep1),7] <- hur_mo_h2.rep1[,4]
rpkm[rownames(hur_mo_h2.rep2),8] <- hur_mo_h2.rep2[,4]
rpkm[rownames(hur_mo_h4.rep1),9] <- hur_mo_h4.rep1[,4]
rpkm[rownames(hur_mo_h4.rep2),10] <- hur_mo_h4.rep2[,4]
rpkm[rownames(hur_mo_h6.rep1),11] <- hur_mo_h6.rep1[,4]
rpkm[rownames(hur_mo_h6.rep2),12] <- hur_mo_h6.rep2[,4]
colnames(rpkm) <- c('ctrl-h2_rep1','ctrl-h2_rep2','ctrl-h4_rep1','ctrl-h4_rep2','ctrl-h6_rep1','ctrl-h6_rep2','hur-h2_rep1','hur-h2_rep2','hur-h4_rep1','hur-h4_rep2','hur-h6_rep1','hur-h6_rep2')


gene.keep <- rownames(rpkm[rpkm[,'ctrl-h2_rep1']+rpkm[,'ctrl-h2_rep2'] > 2,]) ##############要求两小时的时候RPKM>1
gene.keep2 <- rownames(rpkm[rpkm[,'ctrl-h2_rep1']+rpkm[,'ctrl-h2_rep2'] < 2,]) ##############要求两小时的时候RPKM>1
gene.keep3 <- rownames(rpkm[rpkm[,'ctrl-h6_rep1']+rpkm[,'ctrl-h6_rep2'] > 2,]) ##############要求两小时的时候RPKM>1

output <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNseq_ctrl_h26/output_score.txt',sep='\t',stringsAsFactors=F,header=T,row.names=1)
output <- output[complete.cases(output),]

output <- output[intersect(rownames(output),gene.keep),]
fc <- log2(1.5)
fdr <- 0.05

decay <- output[output[,'log2.Fold_change..normalized'] < -log2(1.2) & output[,'q.value.Benjamini.et.al..1995.']<0.05,]

stable <- output[output[,'log2.Fold_change..normalized'] > -log2(1.2)  & output[,'q.value.Benjamini.et.al..1995.']<0.05,]


#################h6/h4 delay的过程
ctrl_h46.out <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq_ctrl_h46/output_score.txt',sep='\t',stringsAsFactors=F,header=T,row.names=1)
mo_h46.out <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-hur-mo-h46/output_score.txt',sep='\t',stringsAsFactors=F,header=T,row.names=1)


ctrl_h46.out <- ctrl_h46.out[complete.cases(ctrl_h46.out),]
mo_h46.out <- mo_h46.out[complete.cases(mo_h46.out),]
####################h4 vs h6 delay

# ctrl_overlay_stable <- intersect(rownames(ctrl_h46.out),rownames(stable))
# mo_overlay_stable <- intersect(rownames(mo_h46.out),rownames(stable))

# ctrl_overlay_dacay <- intersect(rownames(ctrl_h46.out),rownames(decay))
# mo_overlay_dacay <- intersect(rownames(mo_h46.out),rownames(decay))


# x1 <- ecdf(ctrl_h46.out[ctrl_overlay_stable,'log2.Fold_change..normalized'] + seq(0.00000001,by=0.000000001,length.out=length(ctrl_h46.out[ctrl_overlay_stable,'log2.Fold_change..normalized'] )))
# x2 <- ecdf(mo_h46.out[mo_overlay_stable,'log2.Fold_change..normalized'] + seq(0.00000001, by=0.000000001,length.out=length(mo_h46.out[mo_overlay_stable,'log2.Fold_change..normalized'])))


# x3 <- ecdf(ctrl_h46.out[ctrl_overlay_dacay,'log2.Fold_change..normalized'] + seq(0.00000001,by=0.000000001,length.out=length(ctrl_h46.out[ctrl_overlay_dacay,'log2.Fold_change..normalized'] )))
# x4 <- ecdf(mo_h46.out[mo_overlay_dacay,'log2.Fold_change..normalized'] + seq(0.00000001, by=0.000000001,length.out=length(mo_h46.out[mo_overlay_dacay,'log2.Fold_change..normalized'])))

# plot(x1,lwd=1,col="black",main="",xlab="",ylab="", yaxt="n",xlim=c(-10,8), xaxt="n", pch=19, cex=1)
# lines(x2, col="darkblue",do.points=F)
# lines(x3, col="gold",do.points=F)
# lines(x4, col="red",do.points=F)
# abline(v=0)
# dev.off()


# png(paste("h46_delay-cv.png",sep=""),width=600,height=500)
pdf(paste("h46_delay-cv.pdf",sep=""))
par(mgp=c(4,1.5,0), mar=c(7,7,2,2),cex.lab=1.5,cex.axis=1.5)

x1 <- ecdf(ctrl_h46.out[intersect(gene.keep,rownames(ctrl_h46.out)),'log2.Fold_change..normalized'] + seq(0.00000001,by=0.000000001,length.out=length(ctrl_h46.out[intersect(gene.keep,rownames(ctrl_h46.out)),'log2.Fold_change..normalized'] )))

x2 <- ecdf(mo_h46.out[intersect(gene.keep,rownames(mo_h46.out)),'log2.Fold_change..normalized'] + seq(0.00000001,by=0.000000001,length.out=length(mo_h46.out[intersect(gene.keep,rownames(mo_h46.out)),'log2.Fold_change..normalized'] )))

plot(x1,lwd=3,col="black",main="",yaxt="n",xlim=c(-6,6), xaxt="n", pch=19, cex=16, xlab='log2(6h/4h)', ylab='Culmulative frequency')
lines(x2, col="darkgreen",do.points=F,lwd=3)
axis(1,at=seq(-6,6,2),labels=seq(-6,6,2), tcl=-1.3, lwd=3 )
axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2), tcl=-1.3, lwd=3,las=1)
box(lwd=3)
dev.off()








hur_mo_h6.out <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq-hur-mo-h6/output_score.txt',sep='\t',stringsAsFactors=F,header=T,row.names=1)
hur_mo_h6.out <- hur_mo_h6.out[complete.cases(hur_mo_h6.out),]
hur_mo_h6.out <- hur_mo_h6.out[intersect(rownames(hur_mo_h6.out),gene.keep),]

# png(paste("h6_mo-cv.png",sep=""),width=600,height=500)
pdf(paste("h6_mo-cv.pdf",sep=""))
par(mgp=c(4,1.5,0), mar=c(7,7,2,2),cex.lab=1.5,cex.axis=1.5)
hur_h6_overlay_dacay <- intersect(rownames(hur_mo_h6.out),rownames(decay))
hur_h6_overlay_stable <- intersect(rownames(hur_mo_h6.out),rownames(stable))

x1 <- ecdf(hur_mo_h6.out[hur_h6_overlay_dacay,'log2.Fold_change..normalized'] + seq(0.00000001,by=0.000000001,length.out=length(hur_mo_h6.out[hur_h6_overlay_dacay,'log2.Fold_change..normalized'] )))
x2 <- ecdf(hur_mo_h6.out[hur_h6_overlay_stable,'log2.Fold_change..normalized'] + seq(0.00000001,by=0.000000001,length.out=length(hur_mo_h6.out[hur_h6_overlay_stable,'log2.Fold_change..normalized'] )))

plot(x1,lwd=3,col="darkblue",main="",yaxt="n",xlim=c(-2,2), xaxt="n", pch=19, cex=16, xlab='log2(hur_mo/ctrl)', ylab='Culmulative frequency')
lines(x2, col="brown",do.points=F,lwd=3)
axis(1,at=c(-2,-1,0,1,2),labels=c(-2,-1,0,1,2), tcl=-1.3, lwd=3 )
axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2), tcl=-1.3, lwd=3,las=1)
box(lwd=3)
dev.off()