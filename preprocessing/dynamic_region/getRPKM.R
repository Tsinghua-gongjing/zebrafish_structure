setwd("/pnas/yangyg_group/zhangting/SBY_structure-m5C/00-DEGs/HeChuan-group/RPKM")
me <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/00_Combined/0_analysis/hechuan-group/me-all-trans.txt",stringsAsFactors=F)[,2]
ml <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/00_Combined/0_analysis/hechuan-group/ml-all-trans.txt",stringsAsFactors=F)[,2]
s1 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/00_Combined/0_analysis/hechuan-group/s1-all-trans.txt",stringsAsFactors=F)[,2]
s2 <- read.table("//pnas/yangyg_group/zhangting/SBY_structure-m5C/00_Combined/0_analysis/hechuan-group/s2-all-trans.txt",stringsAsFactors=F)[,2]
ze <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/00_Combined/0_analysis/hechuan-group/ze-all-trans.txt",stringsAsFactors=F)[,2]
zl <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/00_Combined/0_analysis/hechuan-group/zl-all-trans.txt",stringsAsFactors=F)[,2]
allgene <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/00_Combined/0_analysis/hechuan-group/total-gene.txt",stringsAsFactors=F)[,2]
# single <- read.table("/pnas/yangyg_group/zhangting/Reference/Zebrafish/ENSEMBLE_79/GeneList/single-transcript-info",stringsAsFactors=F)[,c(6,7)]
# index <- !duplicated(single[,1])
# type <- single[index,]
# rownames(type) <- type[,1]
# table(type[me, 2])

x1 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-8/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x2 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-9/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x3 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-10/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x4 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-11/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x5 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-12/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x6 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-13/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x7 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-14/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x8 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-15/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x9 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-16/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x10 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-20/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x11 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-17/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x12 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-21/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x13 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-18/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x14 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-22/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x15 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-19/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x16 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-23/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x17 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-24/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x18 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-25/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x19 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-26/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x20 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-27/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x21 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-28/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x22 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-29/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x23 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-30/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x24 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-31/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x25 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-32/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x26 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-33/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x27 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-36/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x28 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-37/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x29 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-38/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x30 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-39/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)
x31 <- read.table("/pnas/yangyg_group/zhangting/SBY_structure-m5C/shi-zf-34/analysis/rpkm.txt",row.names=1,header=F,stringsAsFactors=F)


gene.all <- unique(c(rownames(x1),rownames(x2),rownames(x3),rownames(x4),rownames(x5),rownames(x6),rownames(x7),rownames(x8),
rownames(x9),rownames(x10),rownames(x11),rownames(x12),rownames(x13),rownames(x14),rownames(x15),rownames(x16),rownames(x17),
rownames(x18),rownames(x19),rownames(x20),rownames(x21),rownames(x22),rownames(x23),rownames(x24),rownames(x25),rownames(x26),
rownames(x27),rownames(x28),rownames(x29),rownames(x30),rownames(x31),me, ml, s1, s2, ze, zl,allgene))
res <- matrix(0,nrow=length(gene.all),ncol=31)
rownames(res) <- gene.all

############me
res[rownames(x1),1] <- x1[,4]
res[rownames(x2),2] <- x2[,4]
res[rownames(x3),3] <- x3[,4]
res[rownames(x4),4] <- x4[,4]
res[rownames(x5),5] <- x5[,4]
res[rownames(x6),6] <- x6[,4]
res[rownames(x7),7] <- x7[,4]
res[rownames(x8),8] <- x8[,4]
res[rownames(x9),9] <- x9[,4]
res[rownames(x10),10] <- x10[,4]
res[rownames(x11),11] <- x11[,4]
res[rownames(x12),12] <- x12[,4]
res[rownames(x13),13] <- x13[,4]
res[rownames(x14),14] <- x14[,4]
res[rownames(x15),15] <- x15[,4]
res[rownames(x16),16] <- x16[,4]
res[rownames(x17),17] <- x17[,4]
res[rownames(x18),18] <- x18[,4]
res[rownames(x19),19] <- x19[,4]
res[rownames(x20),20] <- x20[,4]
res[rownames(x21),21] <- x21[,4]
res[rownames(x22),22] <- x22[,4]
res[rownames(x23),23] <- x23[,4]
res[rownames(x24),24] <- x24[,4]
res[rownames(x25),25] <- x25[,4]
res[rownames(x26),26] <- x26[,4]
res[rownames(x27),27] <- x27[,4]
res[rownames(x28),28] <- x28[,4]
res[rownames(x29),29] <- x29[,4]
res[rownames(x30),30] <- x30[,4]
res[rownames(x31),31] <- x31[,4]

index <- apply(res, 1, function(fn) all(fn>1))
res <- res[index,]

colnames(res) <- c("egg-DMSO-rep1","egg-NAI-rep1","cell1-DMSO-rep1","cell1-NAI-rep1","cell4-DMSO-rep1","cell4-NAI-rep1",
"cell64-DMSO-rep1","cell64-NAI-rep1","k1-DMSO-rep1","k1-NAI-rep1","sphere-DMSO-rep1","sphere-NAI-rep1","shield-DMSO-rep1",
"shield-NAI-rep1","epiboly-DMSO-rep1","epiboly-NAI-rep1","egg-DMSO-rep2","egg-NAI-rep2","cell1-DMSO-rep2","cell1-NAI-rep2",
"cell4-DMSO-rep2","cell4-NAI-rep2","cell64-DMSO-rep2","cell64-NAI-rep2","k1-DMSO-rep2","k1-NAI-rep2","shield-DMSO-rep2",
"shield-NAI-rep2","epiboly-DMSO-rep2","epiboly-NAI-rep2","sphere-DMSO-rep2")



res.dmso.rpkm <- res[,c('egg-DMSO-rep1','egg-DMSO-rep2','cell1-DMSO-rep1','cell1-DMSO-rep2',
'cell4-DMSO-rep1','cell4-DMSO-rep2','cell64-DMSO-rep1','cell64-DMSO-rep2',
'k1-DMSO-rep1','k1-DMSO-rep2','sphere-DMSO-rep1','sphere-DMSO-rep2','shield-DMSO-rep1',
'shield-DMSO-rep2','epiboly-DMSO-rep1','epiboly-DMSO-rep2')]


res.dmso.com <- cbind(res.dmso.rpkm[,1]+res.dmso.rpkm[,2],res.dmso.rpkm[,3]+res.dmso.rpkm[,4],res.dmso.rpkm[,5]+res.dmso.rpkm[,6],
res.dmso.rpkm[,7]+res.dmso.rpkm[,8],res.dmso.rpkm[,9]+res.dmso.rpkm[,10],res.dmso.rpkm[,11]+res.dmso.rpkm[,12],
res.dmso.rpkm[,13]+res.dmso.rpkm[,14],res.dmso.rpkm[,15]+res.dmso.rpkm[,16])
colnames(res.dmso.com) <- c("egg","cell1","cell4","cell64","k1","sphere","shield","epiboly")

phase1.diff <- res.dmso.com[,"cell64"]/res.dmso.com[,"k1"]
phase2.diff <- res.dmso.com[,"k1"]/res.dmso.com[,"sphere"]
phase3.diff <- res.dmso.com[,"sphere"]/res.dmso.com[,"shield"]

index1 <- phase1.diff > 1.2
index2 <- phase2.diff > 1.2
index3 <- phase3.diff > 1.2

decay.index <- index1 & index2 & index3
decay.gene <- res.dmso.com[decay.index,]
decay.gene <- decay.gene[complete.cases(decay.gene),]

write.table(intersect(rownames(decay.gene),c(me,ml)),"maternal.gold",quote=F,row.names=F,col.names=F)





























# me.rpkm <- res[me,c('egg-DMSO-rep1','egg-DMSO-rep2','cell1-DMSO-rep1','cell1-DMSO-rep2',
# 'cell4-DMSO-rep1','cell4-DMSO-rep2','cell64-DMSO-rep1','cell64-DMSO-rep2',
# 'k1-DMSO-rep1','k1-DMSO-rep2','sphere-DMSO-rep1','sphere-DMSO-rep2','shield-DMSO-rep1',
# 'shield-DMSO-rep2','epiboly-DMSO-rep1','epiboly-DMSO-rep2')]

# me.rpkm.com <- cbind(me.rpkm[,1]+me.rpkm[,2],me.rpkm[,3]+me.rpkm[,4],me.rpkm[,5]+me.rpkm[,6],
# me.rpkm[,7]+me.rpkm[,8],me.rpkm[,9]+me.rpkm[,10],me.rpkm[,11]+me.rpkm[,12],me.rpkm[,13]+me.rpkm[,14],me.rpkm[,15]+me.rpkm[,16])

# colnames(me.rpkm.com) <- c("egg","cell1","cell4","cell64",'1k','sphere','shield','epiboly')

# table(me.rpkm.com[,"egg"]/me.rpkm.com[,"cell64"] >= 1.5)
# table(me.rpkm.com[,"cell64"]/me.rpkm.com[,"sphere"] >= 1.5)

# write.table(me.rpkm.com, "me.rpkm.com",sep="\t",quote=F)



# ml.rpkm <- res[ml,c('egg-DMSO-rep1','egg-DMSO-rep2','cell1-DMSO-rep1','cell1-DMSO-rep2',
# 'cell4-DMSO-rep1','cell4-DMSO-rep2','cell64-DMSO-rep1','cell64-DMSO-rep2',
# 'k1-DMSO-rep1','k1-DMSO-rep2','sphere-DMSO-rep1','sphere-DMSO-rep2','shield-DMSO-rep1',
# 'shield-DMSO-rep2','epiboly-DMSO-rep1','epiboly-DMSO-rep2')]

# ml.rpkm.com <- cbind(ml.rpkm[,1]+ml.rpkm[,2],ml.rpkm[,3]+ml.rpkm[,4],ml.rpkm[,5]+ml.rpkm[,6],
# ml.rpkm[,7]+ml.rpkm[,8],ml.rpkm[,9]+ml.rpkm[,10],ml.rpkm[,11]+ml.rpkm[,12],ml.rpkm[,13]+ml.rpkm[,14],ml.rpkm[,15]+ml.rpkm[,16])

# colnames(ml.rpkm.com) <- c("egg","cell1","cell4","cell64",'1k','sphere','shield','epiboly')

# table(ml.rpkm.com[,"egg"]/ml.rpkm.com[,"cell64"] >= 1.5)
# table(ml.rpkm.com[,"cell64"]/ml.rpkm.com[,"sphere"] >= 1.5)

# table(table(me.rpkm[,'egg-DMSO-rep1']/me.rpkm[,'cell64-DMSO-rep1']>=1.5))

# write.table(ml.rpkm.com, "ml.rpkm.com",sep="\t",quote=F)



# s1.rpkm <- res[s1,c('egg-DMSO-rep1','egg-DMSO-rep2','cell1-DMSO-rep1','cell1-DMSO-rep2',
# 'cell4-DMSO-rep1','cell4-DMSO-rep2','cell64-DMSO-rep1','cell64-DMSO-rep2',
# 'k1-DMSO-rep1','k1-DMSO-rep2','sphere-DMSO-rep1','sphere-DMSO-rep2','shield-DMSO-rep1',
# 'shield-DMSO-rep2','epiboly-DMSO-rep1','epiboly-DMSO-rep2')]

# s1.rpkm.com <- cbind(s1.rpkm[,1]+s1.rpkm[,2],s1.rpkm[,3]+s1.rpkm[,4],s1.rpkm[,5]+s1.rpkm[,6],
# s1.rpkm[,7]+s1.rpkm[,8],s1.rpkm[,9]+s1.rpkm[,10],s1.rpkm[,11]+s1.rpkm[,12],s1.rpkm[,13]+s1.rpkm[,14],s1.rpkm[,15]+s1.rpkm[,16])

# colnames(s1.rpkm.com) <- c("egg","cell1","cell4","cell64",'1k','sphere','shield','epiboly')

# table(s1.rpkm.com[,"egg"]/s1.rpkm.com[,"cell64"] >= 1.5)
# table(s1.rpkm.com[,"cell64"]/s1.rpkm.com[,"sphere"] >= 1.5)
# write.table(s1.rpkm.com, "s1.rpkm.com",sep="\t",quote=F)



# s2.rpkm <- res[s2,c('egg-DMSO-rep1','egg-DMSO-rep2','cell1-DMSO-rep1','cell1-DMSO-rep2',
# 'cell4-DMSO-rep1','cell4-DMSO-rep2','cell64-DMSO-rep1','cell64-DMSO-rep2',
# 'k1-DMSO-rep1','k1-DMSO-rep2','sphere-DMSO-rep1','sphere-DMSO-rep2','shield-DMSO-rep1',
# 'shield-DMSO-rep2','epiboly-DMSO-rep1','epiboly-DMSO-rep2')]

# s2.rpkm.com <- cbind(s2.rpkm[,1]+s2.rpkm[,2],s2.rpkm[,3]+s2.rpkm[,4],s2.rpkm[,5]+s2.rpkm[,6],
# s2.rpkm[,7]+s2.rpkm[,8],s2.rpkm[,9]+s2.rpkm[,10],s2.rpkm[,11]+s2.rpkm[,12],s2.rpkm[,13]+s2.rpkm[,14],s2.rpkm[,15]+s2.rpkm[,16])

# colnames(s2.rpkm.com) <- c("egg","cell1","cell4","cell64",'1k','sphere','shield','epiboly')

# table(s2.rpkm.com[,"egg"]/s2.rpkm.com[,"cell64"] >= 1.5)
# table(s2.rpkm.com[,"cell64"]/s2.rpkm.com[,"sphere"] >= 1.5)

# write.table(s2.rpkm.com, "s2.rpkm.com",sep="\t",quote=F)















































































# ml.rpkm <- res[ml,]
# s1.rpkm <- res[s1,]
# s2.rpkm <- res[s2,]
# ze.rpkm <- res[ze,]
# zl.rpkm <- res[zl,]

# me.egg.rep1 <- table(apply(me.rpkm, 1, function(fn){all(fn[1:2]>1)}))[2]
# me.cell1.rep1 <- table(apply(me.rpkm, 1, function(fn){all(fn[3:4]>1)}))[2]
# me.cell4.rep1 <- table(apply(me.rpkm, 1, function(fn){all(fn[5:6]>1)}))[2]
# me.cell64.rep1 <- table(apply(me.rpkm, 1, function(fn){all(fn[7:8]>1)}))[2]
# me.1k.rep1 <- table(apply(me.rpkm, 1, function(fn){all(fn[9:10]>1)}))[2]
# me.sphere.rep1 <- table(apply(me.rpkm, 1, function(fn){all(fn[11:12]>1)}))[2]
# me.shield.rep1 <- table(apply(me.rpkm, 1, function(fn){all(fn[13:14]>1)}))[2]
# me.epibody.rep1 <- table(apply(me.rpkm, 1, function(fn){all(fn[15:16]>1)}))[2]

# me.egg.rep2 <- table(apply(me.rpkm, 1, function(fn){all(fn[17:18]>1)}))[2]
# me.cell1.rep2 <- table(apply(me.rpkm, 1, function(fn){all(fn[19:20]>1)}))[2]
# me.cell4.rep2 <- table(apply(me.rpkm, 1, function(fn){all(fn[21:22]>1)}))[2]
# me.cell64.rep2 <- table(apply(me.rpkm, 1, function(fn){all(fn[23:24]>1)}))[2]
# me.1k.rep2 <- table(apply(me.rpkm, 1, function(fn){all(fn[25:26]>1)}))[2]
# me.shield.rep2 <- table(apply(me.rpkm, 1, function(fn){all(fn[27:28]>1)}))[2]
# me.epibody.rep2 <- table(apply(me.rpkm, 1, function(fn){all(fn[29:30]>1)}))[2]


# ml.egg.rep1 <- table(apply(ml.rpkm, 1, function(fn){all(fn[1:2]>1)}))[2]
# ml.cell1.rep1 <- table(apply(ml.rpkm, 1, function(fn){all(fn[3:4]>1)}))[2]
# ml.cell4.rep1 <- table(apply(ml.rpkm, 1, function(fn){all(fn[5:6]>1)}))[2]
# ml.cell64.rep1 <- table(apply(ml.rpkm, 1, function(fn){all(fn[7:8]>1)}))[2]
# ml.1k.rep1 <- table(apply(ml.rpkm, 1, function(fn){all(fn[9:10]>1)}))[2]
# ml.sphere.rep1 <- table(apply(ml.rpkm, 1, function(fn){all(fn[11:12]>1)}))[2]
# ml.shield.rep1 <- table(apply(ml.rpkm, 1, function(fn){all(fn[13:14]>1)}))[2]
# ml.epibody.rep1 <- table(apply(ml.rpkm, 1, function(fn){all(fn[15:16]>1)}))[2]

# ml.egg.rep2 <- table(apply(ml.rpkm, 1, function(fn){all(fn[17:18]>1)}))[2]
# ml.cell1.rep2 <- table(apply(ml.rpkm, 1, function(fn){all(fn[19:20]>1)}))[2]
# ml.cell4.rep2 <- table(apply(ml.rpkm, 1, function(fn){all(fn[21:22]>1)}))[2]
# ml.cell64.rep2 <- table(apply(ml.rpkm, 1, function(fn){all(fn[23:24]>1)}))[2]
# ml.1k.rep2 <- table(apply(ml.rpkm, 1, function(fn){all(fn[25:26]>1)}))[2]
# ml.shield.rep2 <- table(apply(ml.rpkm, 1, function(fn){all(fn[27:28]>1)}))[2]
# ml.epibody.rep2 <- table(apply(ml.rpkm, 1, function(fn){all(fn[29:30]>1)}))[2]



# s1.egg.rep1 <- table(apply(s1.rpkm, 1, function(fn){all(fn[1:2]>1)}))[2]
# s1.cell1.rep1 <- table(apply(s1.rpkm, 1, function(fn){all(fn[3:4]>1)}))[2]
# s1.cell4.rep1 <- table(apply(s1.rpkm, 1, function(fn){all(fn[5:6]>1)}))[2]
# s1.cell64.rep1 <- table(apply(s1.rpkm, 1, function(fn){all(fn[7:8]>1)}))[2]
# s1.1k.rep1 <- table(apply(s1.rpkm, 1, function(fn){all(fn[9:10]>1)}))[2]
# s1.sphere.rep1 <- table(apply(s1.rpkm, 1, function(fn){all(fn[11:12]>1)}))[2]
# s1.shield.rep1 <- table(apply(s1.rpkm, 1, function(fn){all(fn[13:14]>1)}))[2]
# s1.epibody.rep1 <- table(apply(s1.rpkm, 1, function(fn){all(fn[15:16]>1)}))[2]

# s1.egg.rep2 <- table(apply(s1.rpkm, 1, function(fn){all(fn[17:18]>1)}))[2]
# s1.cell1.rep2 <- table(apply(s1.rpkm, 1, function(fn){all(fn[19:20]>1)}))[2]
# s1.cell4.rep2 <- table(apply(s1.rpkm, 1, function(fn){all(fn[21:22]>1)}))[2]
# s1.cell64.rep2 <- table(apply(s1.rpkm, 1, function(fn){all(fn[23:24]>1)}))[2]
# s1.1k.rep2 <- table(apply(s1.rpkm, 1, function(fn){all(fn[25:26]>1)}))[2]
# s1.shield.rep2 <- table(apply(s1.rpkm, 1, function(fn){all(fn[27:28]>1)}))[2]
# s1.epibody.rep2 <- table(apply(s1.rpkm, 1, function(fn){all(fn[29:30]>1)}))[2]




# s2.egg.rep1 <- table(apply(s2.rpkm, 1, function(fn){all(fn[1:2]>1)}))[2]
# s2.cell1.rep1 <- table(apply(s2.rpkm, 1, function(fn){all(fn[3:4]>1)}))[2]
# s2.cell4.rep1 <- table(apply(s2.rpkm, 1, function(fn){all(fn[5:6]>1)}))[2]
# s2.cell64.rep1 <- table(apply(s2.rpkm, 1, function(fn){all(fn[7:8]>1)}))[2]
# s2.1k.rep1 <- table(apply(s2.rpkm, 1, function(fn){all(fn[9:10]>1)}))[2]
# s2.sphere.rep1 <- table(apply(s2.rpkm, 1, function(fn){all(fn[11:12]>1)}))[2]
# s2.shield.rep1 <- table(apply(s2.rpkm, 1, function(fn){all(fn[13:14]>1)}))[2]
# s2.epibody.rep1 <- table(apply(s2.rpkm, 1, function(fn){all(fn[15:16]>1)}))[2]

# s2.egg.rep2 <- table(apply(s2.rpkm, 1, function(fn){all(fn[17:18]>1)}))[2]
# s2.cell1.rep2 <- table(apply(s2.rpkm, 1, function(fn){all(fn[19:20]>1)}))[2]
# s2.cell4.rep2 <- table(apply(s2.rpkm, 1, function(fn){all(fn[21:22]>1)}))[2]
# s2.cell64.rep2 <- table(apply(s2.rpkm, 1, function(fn){all(fn[23:24]>1)}))[2]
# s2.1k.rep2 <- table(apply(s2.rpkm, 1, function(fn){all(fn[25:26]>1)}))[2]
# s2.shield.rep2 <- table(apply(s2.rpkm, 1, function(fn){all(fn[27:28]>1)}))[2]
# s2.epibody.rep2 <- table(apply(s2.rpkm, 1, function(fn){all(fn[29:30]>1)}))[2]




# ze.egg.rep1 <- table(apply(ze.rpkm, 1, function(fn){all(fn[1:2]>1)}))[2]
# ze.cell1.rep1 <- table(apply(ze.rpkm, 1, function(fn){all(fn[3:4]>1)}))[2]
# ze.cell4.rep1 <- table(apply(ze.rpkm, 1, function(fn){all(fn[5:6]>1)}))[2]
# ze.cell64.rep1 <- table(apply(ze.rpkm, 1, function(fn){all(fn[7:8]>1)}))[2]
# ze.1k.rep1 <- table(apply(ze.rpkm, 1, function(fn){all(fn[9:10]>1)}))[2]
# ze.sphere.rep1 <- table(apply(ze.rpkm, 1, function(fn){all(fn[11:12]>1)}))[2]
# ze.shield.rep1 <- table(apply(ze.rpkm, 1, function(fn){all(fn[13:14]>1)}))[2]
# ze.epibody.rep1 <- table(apply(ze.rpkm, 1, function(fn){all(fn[15:16]>1)}))[2]

# ze.egg.rep2 <- table(apply(ze.rpkm, 1, function(fn){all(fn[17:18]>1)}))[2]
# ze.cell1.rep2 <- table(apply(ze.rpkm, 1, function(fn){all(fn[19:20]>1)}))[2]
# ze.cell4.rep2 <- table(apply(ze.rpkm, 1, function(fn){all(fn[21:22]>1)}))[2]
# ze.cell64.rep2 <- table(apply(ze.rpkm, 1, function(fn){all(fn[23:24]>1)}))[2]
# ze.1k.rep2 <- table(apply(ze.rpkm, 1, function(fn){all(fn[25:26]>1)}))[2]
# ze.shield.rep2 <- table(apply(ze.rpkm, 1, function(fn){all(fn[27:28]>1)}))[2]
# ze.epibody.rep2 <- table(apply(ze.rpkm, 1, function(fn){all(fn[29:30]>1)}))[2]




# zl.egg.rep1 <- table(apply(zl.rpkm, 1, function(fn){all(fn[1:2]>1)}))[2]
# zl.cell1.rep1 <- table(apply(zl.rpkm, 1, function(fn){all(fn[3:4]>1)}))[2]
# zl.cell4.rep1 <- table(apply(zl.rpkm, 1, function(fn){all(fn[5:6]>1)}))[2]
# zl.cell64.rep1 <- table(apply(zl.rpkm, 1, function(fn){all(fn[7:8]>1)}))[2]
# zl.1k.rep1 <- table(apply(zl.rpkm, 1, function(fn){all(fn[9:10]>1)}))[2]
# zl.sphere.rep1 <- table(apply(zl.rpkm, 1, function(fn){all(fn[11:12]>1)}))[2]
# zl.shield.rep1 <- table(apply(zl.rpkm, 1, function(fn){all(fn[13:14]>1)}))[2]
# zl.epibody.rep1 <- table(apply(zl.rpkm, 1, function(fn){all(fn[15:16]>1)}))[2]

# zl.egg.rep2 <- table(apply(zl.rpkm, 1, function(fn){all(fn[17:18]>1)}))[2]
# zl.cell1.rep2 <- table(apply(zl.rpkm, 1, function(fn){all(fn[19:20]>1)}))[2]
# zl.cell4.rep2 <- table(apply(zl.rpkm, 1, function(fn){all(fn[21:22]>1)}))[2]
# zl.cell64.rep2 <- table(apply(zl.rpkm, 1, function(fn){all(fn[23:24]>1)}))[2]
# zl.1k.rep2 <- table(apply(zl.rpkm, 1, function(fn){all(fn[25:26]>1)}))[2]
# zl.shield.rep2 <- table(apply(zl.rpkm, 1, function(fn){all(fn[27:28]>1)}))[2]
# zl.epibody.rep2 <- table(apply(zl.rpkm, 1, function(fn){all(fn[29:30]>1)}))[2]



