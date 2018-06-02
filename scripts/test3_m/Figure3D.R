# setwd('E://Figure-structure-main//Figure3E')
# sphere <- read.table('E://supplementary//RBP-prediction-gj//sphere-res_extend20.txt',header=T,sep='\t',stringsAsFactors=F)
# shield <- read.table('E://supplementary//RBP-prediction-gj//shield-res_extend20.txt',header=T,sep='\t',stringsAsFactors=F)
# hur <- read.table("E://Figure-structure-supp//FigureS3C-D//hur-mo-h6_down_trans.txt",header=T,sep='\t',stringsAsFactors=F)


# sphere <- cbind(sphere,rowMeans(sphere[,13:59],na.rm=T))
# shield <- cbind(shield,rowMeans(shield[,13:59],na.rm=T))

# rownames(sphere) <- paste(sphere[,1],sphere[,2],sphere[,3],sep='_')
# rownames(shield) <- paste(shield[,1],shield[,2],shield[,3],sep='_')

# com <- intersect(rownames(sphere),rownames(shield))
# com.used <- c()
# for(i in 1:nrow(hur)){
# 	com.used <- c(com.used,com[grepl(hur[i,1], com)])
# }
# sphere.used <- sphere[com,]
# shield.used <- shield[com,]



# decay <- read.table('E://supplementary//decay.txt',sep='\t',stringsAsFactors=F)
# stable <- read.table('E://supplementary//stable.txt',sep='\t',stringsAsFactors=F)
# zygotic <- read.table('E://supplementary//zygotic.txt',sep='\t',stringsAsFactors=F)[,1]
# decay <- rownames(decay)
# stable <- rownames(stable)




# ###############折线图


# maternal.mean <- c()
# stable.mean <- c()
# zygotic.mean <- c()

# for(i in 1:length(decay)){
# 	index <- sphere.used[,1]%in%decay[i]
# 	index2 <- shield.used[,1]%in%decay[i]
# 	maternal.mean <- rbind(maternal.mean, shield.used[index2,13:59]-sphere.used[index,13:59])	
# }

# for(i in 1:length(stable)){
# 	index <- sphere.used[,1]%in%stable[i]
# 	index2 <- shield.used[,1]%in%stable[i]
# 	stable.mean <- rbind(stable.mean, shield.used[index2,13:59]-sphere.used[index,13:59])
# }


# for(i in 1:length(zygotic)){
# 	index <- sphere.used[,1]%in%zygotic[i]
# 	index2 <- shield.used[,1]%in%zygotic[i]
# 	#zygotic.mean <- rbind(zygotic.mean, shield.used[index2,13:59]-sphere.used[index,13:59])
# }


# write.table(maternal.mean, 'maternal-mean.txt',sep='\t',quote=F)
# write.table(stable.mean, 'stable-mean.txt',sep='\t',quote=F)
# write.table(zygotic.mean, 'zygotic-mean.txt',sep='\t',quote=F)\


setwd('/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/Figure-structure-main/Figure3E')

maternal.mean = read.table('maternal-mean.h4_vs_h6_RIP_enrich.h4FC1.5.txt', sep='\t', header=TRUE)
stable.mean = read.table('stable-mean.h4_vs_h6_RIP_enrich.h4FC1.5.txt', sep='\t', header=TRUE)

x <- colMeans(maternal.mean, na.rm=T)
y <- colMeans(stable.mean, na.rm=T)


# png("hur-down-minus.png",width=1000,height=500)
pdf("h4_vs_h6_RIP_enrich.h4FC1.5.txt.pdf",width=12)
par(cex.axis=2,mgp=c(3,3,0),lwd=4,cex.lab=2,mar=c(8,9,2,2))
plot(-1,type="l",lwd=5,col="firebrick",bty="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,47),ylim=c(-0.05,0.05));
lines(1:47,x,lwd=5,col="firebrick",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
points(1:47,x,col="firebrick",cex=1)

lines(1:47,y,lwd=5,col="darkblue",bty="n",xlab="",ylab="",xaxt="n",yaxt="n");
points(1:47,y,col="darkblue",cex=1)
abline(v=21,lwd=2, lty=2)
abline(v=27,lwd=2, lty=2)
axis(1,at=c(0,47),labels=c('',''),tcl=-2, las=1,lwd=4, pos=-0.05)
axis(2,at=c(-0.05,0,0.05),labels=c(-0.05,0,0.05),tcl=-2, las=1,lwd=4, pos=0)
abline(h=0,lty=2,lwd=2)
dev.off()








# #################################################考虑specific
# setwd('E://Figure-structure-main//Figure3I')
# sphere <- read.table('E://supplementary//RBP-prediction-gj//sphere-res_extend20.txt',header=T,sep='\t',stringsAsFactors=F)
# shield <- read.table('E://supplementary//RBP-prediction-gj//shield-res_extend20.txt',header=T,sep='\t',stringsAsFactors=F)
# hur <- read.table("E://Figure-structure-supp//FigureS3C-D//hur-mo-h6_down_trans.txt",header=T,sep='\t',stringsAsFactors=F)


# sphere <- cbind(sphere,rowMeans(sphere[,13:59],na.rm=T))
# shield <- cbind(shield,rowMeans(shield[,13:59],na.rm=T))
# rownames(sphere) <- paste(sphere[,1],sphere[,2],sphere[,3],sep='_')
# rownames(shield) <- paste(shield[,1],shield[,2],shield[,3],sep='_')

# sphere.used <- sphere
# shield.used <- shield


# decay <- read.table('E://supplementary//decay.txt',sep='\t',stringsAsFactors=F)
# stable <- read.table('E://supplementary//stable.txt',sep='\t',stringsAsFactors=F)
# zygotic <- read.table('E://supplementary//zygotic.txt',sep='\t',stringsAsFactors=F)[,1]
# decay <- rownames(decay)
# stable <- rownames(stable)


# sphere.maternal <- c()
# shield.maternal <- c()

# sphere.stable <- c()
# shield.stable <- c()

# for(i in 1:length(decay)){
# 	index <- sphere.used[,1]%in%decay[i]
# 	index2 <- shield.used[,1]%in%decay[i]
# 	sphere.maternal <- rbind(sphere.maternal, sphere.used[index,13:59])
# 	shield.maternal <- rbind(shield.maternal, shield.used[index,13:59])
# }

# for(i in 1:length(stable)){
# 	index <- sphere.used[,1]%in%stable[i]
# 	index2 <- shield.used[,1]%in%stable[i]
# 	sphere.stable <- rbind(sphere.stable, sphere.used[index,13:59])
# 	shield.stable <- rbind(shield.stable, shield.used[index,13:59])
# }


# par(mfrow=c(2,1))
# x <- colMeans(sphere.maternal, na.rm=T)
# y <- colMeans(shield.maternal, na.rm=T)
# plot(1:47, x, type='l', col='red', lwd=2)
# lines(1:47, y, col='green',lwd=2)
# abline(h=0)


# w <- colMeans(sphere.stable, na.rm=T)
# z <- colMeans(shield.stable, na.rm=T)
# plot(1:47, w, type='l', col='red',lwd=2)
# lines(1:47, z, col='green',lwd=2)
# abline(h=0)





# ######################boxplot
# maternal.ratio <- c()
# stable.ratio <- c()
# zygotic.ratio <- c()

# for(i in 1:length(decay)){
# 	index <- sphere.used[,1]%in%decay[i] 
# 	index2 <- shield.used[,1]%in%decay[i]
# 	maternal.ratio <- c(maternal.ratio, shield.used[index2,60]/sphere.used[index,60])
	
# }

# for(i in 1:length(stable)){
# 	index <- sphere.used[,1]%in%stable[i]
# 	index2 <- shield.used[,1]%in%stable[i]
# 	stable.ratio <- c(stable.ratio, shield.used[index2,60]/sphere.used[index,60])	
# }


# for(i in 1:length(zygotic)){
# 	index <- sphere.used[,1]%in%zygotic[i]
# 	index2 <- shield.used[,1]%in%zygotic[i]
# 	zygotic.ratio <- c(zygotic.ratio, shield.used[index2,60]/sphere.used[index,60])	
# }


# boxplot(maternal.ratio,stable.ratio,zygotic.ratio,outline=F)
# abline(h=1,lwd=2,lty=2)

# stable.ratio <- stable.ratio[stable.ratio!="NaN" & stable.ratio!="Inf"]
# maternal.ratio <- maternal.ratio[maternal.ratio!="NaN" & maternal.ratio!="Inf"]
# zygotic.ratio <- zygotic.ratio[zygotic.ratio!="NaN" & zygotic.ratio!="Inf"]
















