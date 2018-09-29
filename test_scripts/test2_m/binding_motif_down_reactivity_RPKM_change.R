setwd('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/Figure-structure-main/Figure2G-H')
abs.sphere <- read.table('./abs/HUR_motif2-sphere.out')
abs.shield <- read.table('./abs/HUR_motif2-sshield.out')

x1 <- 0:40
png("hur-abs-minus.png",width=1000,height=500)
par(cex.axis=2,mgp=c(3,3,0),lwd=4,cex.lab=2,mar=c(8,9,2,2))
plot(-1,type="l",lwd=3,col="black",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,40),ylim=c(-0.1,0.1));
lines(x1,abs.shield-abs.sphere,lwd=5,col="black",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
points(x1,abs.shield-abs.sphere,col="black",cex=1)
abline(v=17,lwd=2, lty=2)
abline(v=23,lwd=2, lty=2)
axis(1,at=c(0,20,40),labels=c(-20,"Site",20),tcl=-2, las=1,lwd=4, pos=-0.1)
axis(2,at=c(-0.1,-0.05,0,0.05,0.1),labels=c(-0.1,-0.05,0,0.05,0.1),tcl=-2, las=1,lwd=4, pos=0)
dev.off()


# up.sphere <- read.table('./increasing/HuR_motif2-sphere.out')
# up.shield <- read.table('./increasing/HuR_motif2-sshield.out')
# png("hur-up-minus.png",width=1000,height=500)
# par(cex.axis=2,mgp=c(3,3,0),lwd=4,cex.lab=2,mar=c(8,9,2,2))
# plot(-1,type="l",lwd=3,col="firebrick",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,40),ylim=c(-0.1,0.1));
# lines(x1,up.shield-up.sphere,lwd=5,col="firebrick",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
# points(x1,up.shield-up.sphere,col="firebrick",cex=1)
# abline(v=17,lwd=2, lty=2)
# abline(v=23,lwd=2, lty=2)
# axis(1,at=c(0,20,40),labels=c(-20,"Site",20),tcl=-2, las=1,lwd=4, pos=-0.1)
# axis(2,at=c(-0.1,-0.05,0,0.05,0.1),labels=c(-0.1,-0.05,0.00,0.05,0.1),tcl=-2, las=1,lwd=4, pos=0)
# dev.off()


down.sphere <- read.table('./decreasing/HuR_motif2-sphere.out')
down.shield <- read.table('./decreasing/HuR_motif2-sshield.out')
png("hur-down-minus.png",width=1000,height=500)
par(cex.axis=2,mgp=c(3,3,0),lwd=4,cex.lab=2,mar=c(8,9,2,2))
plot(-1,type="l",lwd=5,col="darkblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,40),ylim=c(-0.1,0.1));
lines(x1,down.shield-down.sphere,lwd=5,col="darkblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
points(x1,down.shield-down.sphere,col="darkblue",cex=1)
abline(v=17,lwd=2, lty=2)
abline(v=23,lwd=2, lty=2)
axis(1,at=c(0,20,40),labels=c(-20,"Site",20),tcl=-2, las=1,lwd=4, pos=-0.1)
axis(2,at=c(-0.1,-0.05,0,0.05,0.1),labels=c(-0.1,-0.05,0.00,0.05,0.1),tcl=-2, las=1,lwd=4, pos=0)
dev.off()



### merge two lines
pdf("hur-abs-down-minus.pdf", width=12)
par(cex.axis=2,mgp=c(3,3,0),lwd=4,cex.lab=2,mar=c(8,9,2,2))

plot(-1,type="l",lwd=5,col="darkblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,40),ylim=c(-0.1,0.1));
lines(x1,down.shield-down.sphere,lwd=5,col="darkblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
points(x1,down.shield-down.sphere,col="darkblue",cex=1)

# plot(-1,type="l",lwd=3,col="black",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,40),ylim=c(-0.1,0.1));
lines(x1,abs.shield-abs.sphere,lwd=5,col="black",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
points(x1,abs.shield-abs.sphere,col="black",cex=1)

abline(v=17,lwd=2, lty=2)
abline(v=23,lwd=2, lty=2)
axis(1,at=c(0,20,40),labels=c(-20,"Site",20),tcl=-2, las=1,lwd=4, pos=-0.1)
axis(2,at=c(-0.1,-0.05,0,0.05,0.1),labels=c(-0.1,-0.05,0.00,0.05,0.1),tcl=-2, las=1,lwd=4, pos=0)
dev.off()






#######################################基因集合的比较
setwd('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/Figure-structure-main/Figure2G-H')
hur.stable.trans <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/stable-sphere-shield-utr3-fimo_utr3.txt',stringsAsFactors=F)[,1]

hur.up.trans <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/up-sphere-shield-utr3-fimo_utr3.txt',stringsAsFactors=F,sep='\t',row.names=1,header=T)
hur.down.trans <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/down-sphere-shield-utr3-fimo_utr3.txt',stringsAsFactors=F,sep='\t',row.names=1,header=T)

hur.up.trans <- unlist(strsplit(hur.up.trans['HUR_motif2',3],"; "))
hur.down.trans <- unlist(strsplit(hur.down.trans['HUR_motif2',3],"; "))



output <- read.table('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/RNAseq_ctrl_h46/output_score.txt',sep='\t',stringsAsFactors=F,header=T,row.names=1)
output <- output[complete.cases(output),]
fc <- output[,'log2.Fold_change..normalized']
names(fc) <- rownames(output)



# png(paste("hur_utr3-cv.png",sep=""),height=500,width=700)
pdf(paste("hur_utr3-cv.pdf",sep=""))
# par(cex.lab=1.5, cex.axis=1.5, mgp=c(3.5,1,0), mar=c(5,6,5,5))
par(cex.axis=2,mgp=c(3,3,0),lwd=4,cex.lab=2,mar=c(8,9,2,2))
up <- fc[hur.up.trans]
down <- fc[hur.down.trans]
stable <- fc[hur.stable.trans]

change <- fc[unique(c(hur.up.trans,hur.down.trans))]

x1 <- ecdf(up + seq(0.00000001,by=0.000000001,length.out=length(up)))
x2 <- ecdf(down+ seq(0.00000001, by=0.000000001,length.out=length(down)))
x3 <- ecdf(change + seq(0.00000001, by=0.000000001,length.out=length(change)))

plot(x1,lwd=1,col="white",main="",xlab="log2(6h/4h)",ylab="Culmulative Frequency", yaxt="n",xlim=c(-2,2), xaxt="n", pch=19, cex=1.5)
axis(1,at=seq(-2,2,1),labels=rep("",5),tcl=-1.2,lwd=3)
axis(2,at=seq(0,1,0.2),labels=rep("",6),tcl=-1.2,lwd=3)
# lines(x1, col="firebrick",do.points=F, lwd=3)
lines(x2, col="darkblue",do.points=F, lwd=3)
lines(x3, col="black",do.points=F, lwd=3)
# box(lwd=3)
dev.off()

wilcox.test(up, down)$p.value
wilcox.test(up, stable)$p.value
wilcox.test(down, change )$p.value

ks.test(down, stable )$p.value

