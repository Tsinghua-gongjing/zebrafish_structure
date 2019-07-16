setwd('/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/DG91_combine_replicate/Layout/RNAseq2')
ctrl_h2.rep1 <- read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/control_h2-rep1/control_h2-rep1-rpkm.txt",sep='\t',row.names=1,header=F)
ctrl_h2.rep2 <- read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/control_h2-rep2/control_h2-rep2-rpkm.txt",sep='\t',row.names=1,header=F)

ctrl_h4.rep1 <- read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/control_h4-rep1/control_h4-rep1-rpkm.txt",sep='\t',stringsAsFactors=F,header=F,row.names=1)
ctrl_h4.rep2 <- read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/control_h4-rep2/control_h4-rep2-rpkm.txt",sep='\t',stringsAsFactors=F,header=F,row.names=1)

ctrl_h6.rep1 <- read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/control_h6-rep1/control_h6-rep1-rpkm.txt",sep='\t',row.names=1,header=F)
ctrl_h6.rep2 <- read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/control_h6-rep2/control_h6-rep2-rpkm.txt",sep='\t',row.names=1,header=F)

hur_mo_h2.rep1 <- read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/hur_MO-h2_rep1/hur_MO-h2_rep1-rpkm.txt",sep='\t',row.names=1,header=F)
hur_mo_h2.rep2 <- read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/hur_MO-h2_rep2/hur_MO-h2_rep2-rpkm.txt",sep='\t',row.names=1,stringsAsFactors=F)

hur_mo_h4.rep1 <-  read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/hur_MO-h4_rep1/hur_MO-h4_rep1-rpkm.txt",sep='\t',stringsAsFactors=F,header=F,row.names=1)
hur_mo_h4.rep2 <-  read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/hur_MO-h4_rep2/hur_MO-h4_rep2-rpkm.txt",sep='\t',stringsAsFactors=F,header=F,row.names=1)

hur_mo_h6.rep1 <- read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/hur-MO-h6_rep1f/hur-MO-h6_rep1f-rpkm_pe.txt",sep='\t',row.names=1,header=F)
hur_mo_h6.rep2 <- read.table("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/RNAseq2/hur-MO-h6_rep2f/hur-MO-h6_rep2f-rpkm_pe.txt",sep='\t',row.names=1,header=F)

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



##############DEGseq
################DEGseq
#######ctrl_h26
CONTROL<-cbind(rownames(count), count[,c('ctrl-h2_rep1','ctrl-h2_rep2')]);
CASE<-cbind(rownames(count), count[,c('ctrl-h6_rep1','ctrl-h6_rep2')]);

library(DEGseq);
outputdir<-file.path("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/DG91_combine_replicate/Layout/RNAseq2/ctrl_h26/");

layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE));
par(mar=c(2, 2, 2,2));
DEGexp(geneExpMatrix1= CASE, geneCol1=1,expCol1=c(2,3), groupLabel1="h6",geneExpMatrix2=CONTROL,geneCol2=1,expCol2=c(2,3),groupLabel2="h2",outputDir=outputdir,method="MARS");


#######ctrl_h46
CONTROL<-cbind(rownames(count), count[,c('ctrl-h4_rep1','ctrl-h4_rep2')]);
CASE<-cbind(rownames(count), count[,c('ctrl-h6_rep1','ctrl-h6_rep2')]);

library(DEGseq);
outputdir<-file.path("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/DG91_combine_replicate/Layout/RNAseq2/ctrl_h46/");

layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE));
par(mar=c(2, 2, 2,2));
DEGexp(geneExpMatrix1= CASE, geneCol1=1,expCol1=c(2,3), groupLabel1="h6",geneExpMatrix2=CONTROL,geneCol2=1,expCol2=c(2,3),groupLabel2="h2",outputDir=outputdir,method="MARS");






#######hur_mo_h26
CONTROL<-cbind(rownames(count), count[,c('hur-h2_rep1','hur-h2_rep2')]);
CASE<-cbind(rownames(count), count[,c('hur-h6_rep1','hur-h6_rep2')]);

library(DEGseq);
outputdir<-file.path("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/DG91_combine_replicate/Layout/RNAseq2/hur_mo_h26/");

layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE));
par(mar=c(2, 2, 2,2));
DEGexp(geneExpMatrix1= CASE, geneCol1=1,expCol1=c(2,3), groupLabel1="h6",geneExpMatrix2=CONTROL,geneCol2=1,expCol2=c(2,3),groupLabel2="h2",outputDir=outputdir,method="MARS");




#######hur_mo_h46
CONTROL<-cbind(rownames(count), count[,c('hur-h4_rep1','hur-h4_rep2')]);
CASE<-cbind(rownames(count), count[,c('hur-h6_rep1','hur-h6_rep2')]);

library(DEGseq);
outputdir<-file.path("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/DG91_combine_replicate/Layout/RNAseq2/hur_mo_h46/");

layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE));
par(mar=c(2, 2, 2,2));
DEGexp(geneExpMatrix1= CASE, geneCol1=1,expCol1=c(2,3), groupLabel1="h6",geneExpMatrix2=CONTROL,geneCol2=1,expCol2=c(2,3),groupLabel2="h2",outputDir=outputdir,method="MARS");



#######hur_mo_h6/ctr_h6
CONTROL<-cbind(rownames(count), count[,c('ctrl-h6_rep1','ctrl-h6_rep2')]);
CASE<-cbind(rownames(count), count[,c('hur-h6_rep1','hur-h6_rep2')]);

library(DEGseq);
outputdir<-file.path("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/DG91_combine_replicate/Layout/RNAseq2/h6/");

layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE));
par(mar=c(2, 2, 2,2));
DEGexp(geneExpMatrix1= CASE, geneCol1=1,expCol1=c(2,3), groupLabel1="h6",geneExpMatrix2=CONTROL,geneCol2=1,expCol2=c(2,3),groupLabel2="h2",outputDir=outputdir,method="MARS");




#######hur_mo_h4/ctr_h4
CONTROL<-cbind(rownames(count), count[,c('ctrl-h4_rep1','ctrl-h4_rep2')]);
CASE<-cbind(rownames(count), count[,c('hur-h4_rep1','hur-h4_rep2')]);

library(DEGseq);
outputdir<-file.path("/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/DG91_combine_replicate/Layout/RNAseq2/h4");

layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE));
par(mar=c(2, 2, 2,2));
DEGexp(geneExpMatrix1= CASE, geneCol1=1,expCol1=c(2,3), groupLabel1="h6",geneExpMatrix2=CONTROL,geneCol2=1,expCol2=c(2,3),groupLabel2="h2",outputDir=outputdir,method="MARS");






###################定义maternal/stable集合
setwd('/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/DG91_combine_replicate/Layout/RNAseq2')
output <- read.table('/pnas/yangyg_group/hanyn/zhangting/zebrafish_paris/DG91_combine_replicate/Layout/RNAseq2/ctrl_h26/output_score.txt',sep='\t',stringsAsFactors=F,header=T,row.names=1)
output <- output[complete.cases(output),]
gene.keep <- rownames(rpkm[rpkm[,'ctrl-h2_rep1']+rpkm[,'ctrl-h2_rep2'] > 2,]) ##############要求两小时的时候RPKM>1
gene.keep2 <- rownames(rpkm[rpkm[,'ctrl-h2_rep1']+rpkm[,'ctrl-h2_rep2'] < 2,]) ##############要求两小时的时候RPKM>1
gene.keep3 <- rownames(rpkm[rpkm[,'ctrl-h6_rep1']+rpkm[,'ctrl-h6_rep2'] > 2,]) ##############要求两小时的时候RPKM>1


write.table(gene.keep,"gene-keep.txt",quote=F,row.names=F,col.names=F)
write.table(intersect(gene.keep3,gene.keep2),"zygotic.txt",quote=F,row.names=F,col.names=F)


output <- output[intersect(rownames(output),gene.keep),]
fc <- log2(1.5)
fdr <- 0.05

decay <- output[output[,'log2.Fold_change..normalized'] < -log2(1.2) & output[,'q.value.Benjamini.et.al..1995.']<0.05,]

stable <- output[output[,'log2.Fold_change..normalized'] > -log2(1.2)  & output[,'q.value.Benjamini.et.al..1995.']<0.05,]



write.table(decay,'decay.txt',sep='\t',quote=F)
write.table(stable,'stable.txt',sep='\t',quote=F)








