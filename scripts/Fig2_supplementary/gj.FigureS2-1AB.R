#######################################
library(data.table)
options(stringsAsFactors=F)
getData <- function(fn){
	tmp.data <- fread(fn,header=F, sep="\t")
	tmp.data <- data.frame(tmp.data)
	rownames(tmp.data) <- tmp.data[,1]
	tmp.data <- tmp.data[,-1]
}

setwd("/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/Figure")
#############################################
total <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/smooth-bigmatrix-stage6.txt"
up <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/matrix.variable-up.txt"
down <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/matrix.variable-down.txt"

up.data <- getData(up)
down.data <- getData(down)
total.data <- getData(total)
res1 <- c() 
res2 <- c()

num <- dim(total.data)[2]-1
for(i in seq(1,num,1)){
	
	stage <- total.data[,i:(i+1)]
	index <- apply(stage, 1, function(fn){all(fn != "NULL")})
	total.num <- dim(stage[index,])[1]
	
	up.index <- !is.na(up.data[,i])
	up.num <- dim(up.data[up.index,])[1]
	down.index <- !is.na(down.data[,i])
	down.num <- dim(down.data[down.index,])[1]
	
	res1 <- cbind( res1, c(total.num-(up.num+down.num),up.num+down.num))
	res2 <- cbind( res2, c(down.num, up.num) )
}
max(res1)
max(res2)
pdf("Figure2A-01_001.pdf", width=12)
#png("Figure2A-01_001.png", height=500)
par(lwd=3, cex.axis=1.5, cex.lab=2, las=2, mar=c(10,5,3,10))	
barplot(res1, col=c("#9BC3E3", "#4462A5"),ylim=c(0, 1000000), xlim=c(0,2), xlab='', ylab='Count (10bp Tiles)', yaxt='n', width=0.3, space=0.4, names.arg=c("Egg->1 Cell","1 Cell->4Cell","4 Cell->64 Cell","64 Cell->Sphere", "Sphere->Shield"),legend.text=c("Stable","Change"), args.legend = list(x=0.8, y=900000, border="white"))
axis(2, at=c(0, 500000 ,1000000), labels=c(0,5,10), las=1,lwd=3, tcl=-1,pos=0 )
dev.off()
pdf("Figure2B-01_001.pdf", width=12)
#png("Figure2B-01_001.png")
par(lwd=3, cex.axis=1.5, cex.lab=2, las=2, mar=c(10,5,3,10))		
barplot(res2, col=c("#74B72E", "#DD3E4B"), lwd=3, ylim=c(0,90000), xlim=c(0,2), xlab='', ylab='Count (10bp Tiles)', yaxt='n', width=0.3, space=0.4, names.arg=c("Egg->1Cell","1 Cell->4Cell","4 Cell->64 Cell","64 Cell->Sphere", "Sphere->Shield"),legend.text=c("Decreasing","Increasing"), args.legend = list("toprigt",border="white"))
axis(2,at=c(0,45000,90000), labels=c(0,4.5,9), pos=0,las=1, lwd=3, tcl=-1 )
dev.off()



# #############################################
# up <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_005/matrix.variable-up.txt"
# down <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_005/matrix.variable-down.txt"
# up.data <- getData(up)
# down.data <- getData(down)
# total.data <- getData(total)
# res1 <- c() 
# res2 <- c()

# num <- dim(total.data)[2]-1
# for(i in seq(1,num,1)){
	
# 	stage <- total.data[,i:(i+1)]
# 	index <- apply(stage, 1, function(fn){all(fn != "NULL")})
# 	total.num <- dim(stage[index,])[1]
	
# 	up.index <- !is.na(up.data[,i])
# 	up.num <- dim(up.data[up.index,])[1]
# 	down.index <- !is.na(down.data[,i])
# 	down.num <- dim(down.data[down.index,])[1]
	
# 	res1 <- cbind( res1, c(total.num-(up.num+down.num),up.num+down.num))
# 	res2 <- cbind( res2, c(down.num, up.num) )
# }
# max(res1)
# max(res2)


# png("Figure2A-01_005.png", height=500)
# par(lwd=3, cex.axis=1.5, cex.lab=2, las=2, mar=c(10,5,3,10))	
# barplot(res1, col=c("#9BC3E3", "#4462A5"),ylim=c(0, 1000000), xlim=c(0,2), xlab='', ylab='Count (10bp Tiles)', yaxt='n', width=0.3, space=0.4, names.arg=c("Egg->1 Cell","1 Cell->4Cell","4 Cell->64 Cell","64 Cell->1 K","1 K->Sphere", "Sphere->Shield"),legend.text=c("Stable","Change"), args.legend = list(x=0.8, y=900000, border="white"))
# axis(2, at=c(0, 500000 ,1000000), labels=c(0,5,10), las=1,lwd=3, tcl=-1,pos=0 )
# dev.off()

# png("Figure2B-01_005.png")
# par(lwd=3, cex.axis=1.5, cex.lab=2, las=2, mar=c(10,5,3,10))		
# barplot(res2, col=c("#74B72E", "#DD3E4B"), lwd=3, ylim=c(0,180000), xlim=c(0,2), xlab='', ylab='Count (10bp Tiles)', yaxt='n', width=0.3, space=0.4, names.arg=c("Egg->1Cell","1 Cell->4Cell","4 Cell->64 Cell","64 Cell->1 K","1 K->Sphere", "Sphere->Shield"),legend.text=c("Decreasing","Increasing"), args.legend = list("toprigt",border="white"))
# axis(2,at=c(0,90000,180000), labels=c(0,9,1.8), pos=0,las=1, lwd=3, tcl=-1 )
# dev.off()


# #############################################
# up <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_001/matrix.variable-up.txt"
# down <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_001/matrix.variable-down.txt"
# up.data <- getData(up)
# down.data <- getData(down)
# total.data <- getData(total)
# res1 <- c() 
# res2 <- c()

# num <- dim(total.data)[2]-1
# for(i in seq(1,num,1)){
	
# 	stage <- total.data[,i:(i+1)]
# 	index <- apply(stage, 1, function(fn){all(fn != "NULL")})
# 	total.num <- dim(stage[index,])[1]
	
# 	up.index <- !is.na(up.data[,i])
# 	up.num <- dim(up.data[up.index,])[1]
# 	down.index <- !is.na(down.data[,i])
# 	down.num <- dim(down.data[down.index,])[1]
	
# 	res1 <- cbind( res1, c(total.num-(up.num+down.num),up.num+down.num))
# 	res2 <- cbind( res2, c(down.num, up.num) )
# }
# max(res1)
# max(res2)

# png("Figure2A-005_001.png", height=500)
# par(lwd=3, cex.axis=1.5, cex.lab=2, las=2, mar=c(10,5,3,10))	
# barplot(res1, col=c("#9BC3E3", "#4462A5"),ylim=c(0, 1000000), xlim=c(0,2), xlab='', ylab='Count (10bp Tiles)', yaxt='n', width=0.3, space=0.4, names.arg=c("Egg->1 Cell","1 Cell->4Cell","4 Cell->64 Cell","64 Cell->1 K","1 K->Sphere", "Sphere->Shield"),legend.text=c("Stable","Change"), args.legend = list(x=0.8, y=900000, border="white"))
# axis(2, at=c(0, 500000 ,1000000), labels=c(0,5,10), las=1,lwd=3, tcl=-1,pos=0 )
# dev.off()

# png("Figure2B-005_001.png")
# par(lwd=3, cex.axis=1.5, cex.lab=2, las=2, mar=c(10,5,3,10))		
# barplot(res2, col=c("#74B72E", "#DD3E4B"), lwd=3, ylim=c(0,120000), xlim=c(0,2), xlab='', ylab='Count (10bp Tiles)', yaxt='n', width=0.3, space=0.4, names.arg=c("Egg->1Cell","1 Cell->4Cell","4 Cell->64 Cell","64 Cell->1 K","1 K->Sphere", "Sphere->Shield"),legend.text=c("Decreasing","Increasing"), args.legend = list("toprigt",border="white"))
# axis(2,at=c(0,60000,120000), labels=c(0,6,12), pos=0,las=1, lwd=3, tcl=-1 )
# dev.off()



# #############################################
# up <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_005/matrix.variable-up.txt"
# down <- "/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_005/matrix.variable-down.txt"
# up.data <- getData(up)
# down.data <- getData(down)
# total.data <- getData(total)
# res1 <- c() 
# res2 <- c()

# num <- dim(total.data)[2]-1
# for(i in seq(1,num,1)){
	
# 	stage <- total.data[,i:(i+1)]
# 	index <- apply(stage, 1, function(fn){all(fn != "NULL")})
# 	total.num <- dim(stage[index,])[1]
	
# 	up.index <- !is.na(up.data[,i])
# 	up.num <- dim(up.data[up.index,])[1]
# 	down.index <- !is.na(down.data[,i])
# 	down.num <- dim(down.data[down.index,])[1]
	
# 	res1 <- cbind( res1, c(total.num-(up.num+down.num),up.num+down.num))
# 	res2 <- cbind( res2, c(down.num, up.num) )
# }
# max(res1)
# max(res2)

# png("Figure2A-005_005.png", height=500)
# par(lwd=3, cex.axis=1.5, cex.lab=2, las=2, mar=c(10,5,3,10))	
# barplot(res1, col=c("#9BC3E3", "#4462A5"),ylim=c(0, 1000000), xlim=c(0,2), xlab='', ylab='Count (10bp Tiles)', yaxt='n', width=0.3, space=0.4, names.arg=c("Egg->1 Cell","1 Cell->4Cell","4 Cell->64 Cell","64 Cell->1 K","1 K->Sphere", "Sphere->Shield"),legend.text=c("Stable","Change"), args.legend = list(x=0.8, y=900000, border="white"))
# axis(2, at=c(0, 500000 ,1000000), labels=c(0,5,10), las=1,lwd=3, tcl=-1,pos=0 )
# dev.off()

# png("Figure2B-005_005.png")
# par(lwd=3, cex.axis=1.5, cex.lab=2, las=2, mar=c(10,5,3,10))		
# barplot(res2, col=c("#74B72E", "#DD3E4B"), lwd=3, ylim=c(0,220000), xlim=c(0,2), xlab='', ylab='Count (10bp Tiles)', yaxt='n', width=0.3, space=0.4, names.arg=c("Egg->1Cell","1 Cell->4Cell","4 Cell->64 Cell","64 Cell->1 K","1 K->Sphere", "Sphere->Shield"),legend.text=c("Decreasing","Increasing"), args.legend = list("toprigt",border="white"))
# axis(2,at=c(0,110000,220000), labels=c(0,11,22), pos=0,las=1, lwd=3, tcl=-1 )
# dev.off()








