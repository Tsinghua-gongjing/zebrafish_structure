setwd('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/Figure-structure-main/Figure3E')

maternal.mean = read.table('maternal-mean.new.h4_vs_h6_RIP_enrich.h4FC1.5.txt', sep='\t', header=TRUE)
stable.mean = read.table('stable-mean.new.h4_vs_h6_RIP_enrich.h4FC1.5.txt', sep='\t', header=TRUE)


x <- colMeans(maternal.mean, na.rm=T)
y <- colMeans(stable.mean, na.rm=T)
Pval <- c()
for(i in 1:dim(maternal.mean)[2]){
	# tsum = t.test(maternal.mean[,i], stable.mean[,i])
	tsum = wilcox.test(maternal.mean[,i], stable.mean[,i])
	Pval <- c(Pval, tsum$p.value)
}
print(Pval)


# png("hur-down-minus.png",width=1000,height=500)
# pdf("h4_vs_h6_RIP_enrich.h4FC1.5.txt.new.pdf",width=12)
# pdf("hur-down-minus.pdf",width=12)
pdf("/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/FS4.h4_vs_h6.iCLIPbind.reactivity_meta.pdf",width=12)
par(cex.axis=2,mgp=c(3,3,0),lwd=4,cex.lab=2,mar=c(8,9,2,2))
plot(-1,type="l",lwd=5,col="firebrick",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,47),ylim=c(-0.06,0.06));
lines(1:47,x,lwd=5,col="firebrick",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
points(1:47,x,col="firebrick",cex=1)

lines(1:47,y,lwd=5,col="darkblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
points(1:47,y,col="darkblue",cex=1)
abline(v=21,lwd=2, lty=2)
abline(v=27,lwd=2, lty=2)
axis(1,at=c(0,47),labels=c('',''),tcl=-2, las=1,lwd=4, pos=-0.05)
axis(2,at=c(-0.06,0,0.06),labels=c(-0.06,0,0.06),tcl=-2, las=1,lwd=4, pos=0)
abline(h=0,lty=2,lwd=2)

print(dim(maternal.mean))
print(dim(stable.mean))

sig_pos = c()
sig_val = c()
for(i in 1:47){
	if(Pval[i] <= 0.05){
		sig_pos <- c(sig_pos, i)
		sig_val <- c(sig_val, Pval[i])
		if(Pval[i] <= 0.01){
			text(i, y[i]+0.01, '**', cex=2, srt = 90)
		}else{
				text(i, y[i]+0.01, '*', cex=2, srt = 90)
			}
		
	}
}
print(sig_pos)
print(sig_val)


dev.off()