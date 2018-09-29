
getData <- function(fn, element){
	
	dataSet <- fread(fn, header=F,sep="\t")
	dataSet <- data.frame(dataSet)
	rowNames <- paste(dataSet[,1],dataSet[,2],dataSet[,3],sep="_")
	dataSet <- dataSet[,-c(1:3)]
	dataSet <- data.frame(dataSet)
	#rownames(dataSet) <- rowNames
	dataSet
}

##服务器的地址
setwd("/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/Figure")
x1 <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/Figure/smooth/egg-smooth_anno.bed"
x2 <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/Figure/smooth/cell1-smooth_anno.bed"
x3 <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/Figure/smooth/cell4-smooth_anno.bed"
x4 <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/Figure/smooth/cell64-smooth_anno.bed"
x5 <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/Figure/smooth/sphere-smooth_anno.bed"
x6 <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/Figure/smooth/shield-smooth_anno.bed"

library(data.table)
options(stringsAsFactors=F)
plotData <- list()
plotData[[1]] <- getData(x1)
plotData[[2]] <- getData(x2)
plotData[[3]] <- getData(x3)
plotData[[4]] <- getData(x4)
plotData[[5]] <- getData(x5)
plotData[[6]] <- getData(x6)


#####high median low

drawFreq <- function(dataSet, picName){
	y <- length(dataSet)
	#png("utr3_freq.png", width=2800, height=1500)
	pdf("utr3_freq.pdf", width=10)
	par(mfrow=c(2,y/2),cex.axis=5, cex.lab=8, lwd=5, las=1, mar=c(20,15,6,6), mgp=c(15,6,0))
	for(i in 1:y){
		index <- grepl("utr3", dataSet[[i]][,6])
		tmp.res <- apply(dataSet[[i]][index,], 1, function(fn){ unlist(strsplit(fn,","))[1]})
		trans <- sapply(strsplit(rownames(dataSet[[i]][index, ]),'_'),function(fn){ paste(fn[1],fn[2],sep="_")})
		transLen <- length(unique(trans))
		hist(as.numeric(tmp.res), xlim=c(0,1),ylim=c(0,50000), main='',xlab=picName[i], ylab="", xaxt='n',yaxt='n',breaks=20, col="#3B3751")
		axis(1, at=seq(0,1,0.2), labels=seq(0,1,0.2), tcl=-3, lwd=7, pos=0)
		axis(2, at=c(0,25000,50000), labels=c(0,2.5,5), tcl=-3,las=1, lwd=7, pos=0)
		legend('topright',c(paste("W=",length(tmp.res)), paste('T=',transLen)), bty="n", cex=5)
	}
	dev.off()
	
	#png("cds_freq.png", width=2800, height=1500)
	pdf("cds_freq.pdf", width=10)
	par(mfrow=c(2,y/2),cex.axis=5, cex.lab=8, lwd=5, las=1, mar=c(20,15,6,6), mgp=c(15,6,0))
	for(i in 1:y){
		index <- grepl("cds", dataSet[[i]][,6])
		tmp.res <- apply(dataSet[[i]][index,], 1, function(fn){ unlist(strsplit(fn,","))[1]})
		trans <- sapply(strsplit(rownames(dataSet[[i]][index, ]),'_'),function(fn){ paste(fn[1],fn[2],sep="_")})
		transLen <- length(unique(trans))
		hist(as.numeric(tmp.res), xlim=c(0,1),ylim=c(0,200000), main='',xlab=picName[i], ylab="", xaxt='n',yaxt='n',breaks=20, col="#3B3751")
		
		axis(1, at=seq(0,1,0.2), labels=seq(0,1,0.2), tcl=-3, lwd=7, pos=0)
		axis(2, at=c(0,100000,200000), labels=c(0,10,20), tcl=-3,las=1, lwd=7, pos=0)
		legend('topright',c(paste("W=",length(tmp.res)), paste('T=',transLen)), bty="n", cex=5)
	}
	dev.off()
	
	#png("utr5_freq.png", width=2800, height=1500)
	pdf("utr5_freq.pdf", width=10)
	par(mfrow=c(2,y/2),cex.axis=5, cex.lab=8, lwd=5, las=1, mar=c(20,15,6,6), mgp=c(15,6,0))
	for(i in 1:y){
		index <- grepl("utr5", dataSet[[i]][,6])
		tmp.res <- apply(dataSet[[i]][index,], 1, function(fn){ unlist(strsplit(fn,","))[1]})
		trans <- sapply(strsplit(rownames(dataSet[[i]][index, ]),'_'),function(fn){ paste(fn[1],fn[2],sep="_")})
		transLen <- length(unique(trans))
		hist(as.numeric(tmp.res), xlim=c(0,1),ylim=c(0,30000), main='',xlab=picName[i], ylab="", xaxt='n',yaxt='n',breaks=20, col="#3B3751")
		
		axis(1, at=seq(0,1,0.2), labels=seq(0,1,0.2), tcl=-3, lwd=7, pos=0)
		axis(2, at=c(0,15000,30000), labels=c(0,1.5,3), tcl=-3,las=1, lwd=7, pos=0)
		legend('topright',c(paste("W=",length(tmp.res)), paste('T=',transLen)), bty="n", cex=5)
	}
	dev.off()
	
	
	
	
}

drawFreq(plotData, c("Egg", "1 Cell","4 Cell","64 Cell", "Sphere", "Shield") )

drawHML <- function(dataSet, high,low){
	h <- high
	l <- low
	y <- length(dataSet)
	res <- c()
	for(i in 1:y){	
		tmp.res <- apply(data.frame(dataSet[[i]][,1]), 1 , function(fn){ unlist(strsplit(fn,","))[1]})
		tmp <- as.numeric(tmp.res)
		t.h.n <- length(tmp[tmp>=h])/length(tmp)
		t.l.n <- length(tmp[tmp<=l])/length(tmp)
		t.m.n <- 1-t.h.n-t.l.n
		res <- cbind(res,c(t.l.n,t.m.n,t.h.n))
	}
	#png("transcript-bar.png")
	pdf("transcript-bar.pdf")
	par(lwd=3, cex.axis=1.5, cex.lab=2.5, mar=c(5,6,6,6), las=2)	
	barplot(res, col=c("navy", "limegreen", "firebrick2"), lwd=3, ylim=c(0,1), xlim=c(0,2), xlab='', ylab='%Total', yaxt='n', width=0.25, names.arg=c("Egg", "1 Cell","4 Cell","64 Cell","Sphere", "Shield"),legend.text=c("Low","Median","High"), args.legend = list(x=2, y=0.7,bg="white"))
	axis(2,at=seq(0,1,0.2), labels=seq(0,1,0.2),pos=0,las=1, lwd=3, tcl=-1 )
	dev.off()

	res <- c()
	for(i in 1:y){	
		index <- grepl("utr3", dataSet[[i]][,6])
		tmp.res <- apply(data.frame(dataSet[[i]][index,1]), 1 , function(fn){ unlist(strsplit(fn,","))[1]})
		tmp <- as.numeric(tmp.res)
		t.h.n <- length(tmp[tmp>=h])/length(tmp)
		t.l.n <- length(tmp[tmp<=l])/length(tmp)
		t.m.n <- 1-t.h.n-t.l.n
		res <- cbind(res,c(t.l.n,t.m.n,t.h.n))
	}
	#png("utr3-bar.png")
	pdf("utr3-bar.pdf")
	par(lwd=3, cex.axis=1.5, cex.lab=2.5, mar=c(5,6,6,6), las=2)	
	barplot(res, col=c("navy", "limegreen", "firebrick2"), lwd=3, ylim=c(0,1), xlim=c(0,2), xlab='', ylab='%Total', yaxt='n', width=0.25, names.arg=c("Egg", "1 Cell","4 Cell","64 Cell","Sphere", "Shield"),legend.text=c("Low","Median","High"), args.legend = list(x=2, y=0.7,bg="white"))
	axis(2,at=seq(0,1,0.2), labels=seq(0,1,0.2),pos=0,las=1, lwd=3, tcl=-1 )
	dev.off()
	
	res <- c()
	for(i in 1:y){	
		index <- grepl("utr5", dataSet[[i]][,6])
		tmp.res <- apply(data.frame(dataSet[[i]][index,1]), 1 , function(fn){ unlist(strsplit(fn,","))[1]})
		tmp <- as.numeric(tmp.res)
		t.h.n <- length(tmp[tmp>=h])/length(tmp)
		t.l.n <- length(tmp[tmp<=l])/length(tmp)
		t.m.n <- 1-t.h.n-t.l.n
		res <- cbind(res,c(t.l.n,t.m.n,t.h.n))
	}
	#png("utr5-bar.png")
	pdf("utr5-bar.pdf")
	par(lwd=3, cex.axis=1.5, cex.lab=2.5, mar=c(5,6,6,6), las=2)	
	barplot(res, col=c("navy", "limegreen", "firebrick2"), lwd=3, ylim=c(0,1), xlim=c(0,2), xlab='', ylab='%Total', yaxt='n', width=0.25, names.arg=c("Egg", "1 Cell","4 Cell","64 Cell","Sphere", "Shield"),legend.text=c("Low","Median","High"), args.legend = list(x=2, y=0.7,bg="white"))
	axis(2,at=seq(0,1,0.2), labels=seq(0,1,0.2),pos=0,las=1, lwd=3, tcl=-1 )
	dev.off()
	
	
	res <- c()
	for(i in 1:y){	
		index <- grepl("cds", dataSet[[i]][,6])
		tmp.res <- apply(data.frame(dataSet[[i]][index,1]), 1 , function(fn){ unlist(strsplit(fn,","))[1]})
		tmp <- as.numeric(tmp.res)
		t.h.n <- length(tmp[tmp>=h])/length(tmp)
		t.l.n <- length(tmp[tmp<=l])/length(tmp)
		t.m.n <- 1-t.h.n-t.l.n
		res <- cbind(res,c(t.l.n,t.m.n,t.h.n))
	}
	#png("cds-bar.png")
	pdf("cds-bar.pdf")
	par(lwd=3, cex.axis=1.5, cex.lab=2.5, mar=c(5,6,6,6), las=2)	
	barplot(res, col=c("navy", "limegreen", "firebrick2"), lwd=3, ylim=c(0,1), xlim=c(0,2), xlab='', ylab='%Total', yaxt='n', width=0.25, names.arg=c("Egg", "1 Cell","4 Cell","64 Cell","Sphere", "Shield"),legend.text=c("Low","Median","High"), args.legend = list(x=2, y=0.7,bg="white"))
	axis(2,at=seq(0,1,0.2), labels=seq(0,1,0.2),pos=0,las=1, lwd=3, tcl=-1 )
	dev.off()	
}
drawHML(plotData,0.4, 0.2)


















