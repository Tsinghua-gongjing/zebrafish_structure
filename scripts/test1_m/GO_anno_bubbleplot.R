library(ggplot2)
data <- read.table("bottom-bp_used.txt",sep='\t',header=T,stringsAsFactors = F)

data[,1] <- matrix(unlist(strsplit(data[,1],"~")),byrow=T,nrow(data))[,2]
data <- cbind(data, 1:nrow(data))
#data
#data2 <- data
colnames(data) <- ""
# plot.data <- rbind(data[,c(1,2,3,ncol(data))],data[,c(1,4,5,ncol(data))],data[,c(1,6,7,ncol(data))],data[,c(1,8,9,ncol(data))],data[,c(1,10,11,ncol(data))],data[,c(1,12,13,ncol(data))],data[,c(1,14,15,ncol(data))])
#plot.data <- rbind(data[,c(1,2,3,ncol(data))],data[,c(1,4,5,ncol(data))],data[,c(1,6,7,ncol(data))],data[,c(1,8,9,ncol(data))],data[,c(1,10,11,ncol(data))],data[,c(1,12,13,ncol(data))])
plot.data <- rbind(data[,c(1,2,3,4,ncol(data))],data[,c(1,5,6,7,ncol(data))],data[,c(1,8,9,10,ncol(data))],data[,c(1,11,12,13,ncol(data))],data[,c(1,14,15,16,ncol(data))],data[,c(1,17,18,19,ncol(data))])
# plot.data$stage <- rep(c("0-Egg","1-Cell","4-cells","64-cells","K1","Sphere","Shield"),each=nrow(data))
plot.data$stage <- rep(c("Egg","1-Cell","4-Cell","64-Cell","Sphere","Shield"),each=nrow(data))
colnames(plot.data) <- c("term","count","p","fdr","termorder","stage")
#plot.data$p <- -log10(plot.data$p)
#plot.data$fdr <- -log10(plot.data$fdr)

library(Hmisc)
plot.data$term <-  capitalize(plot.data$term)
plot.data$term <- factor(plot.data$term,levels=plot.data[order(plot.data$termorder, decreasing = T),1])
#plot.data$term <- factor(plot.data$term,levels=plot.data[order(plot.data$fdr, decreasing = T),1])
plot.data$stage <- factor(plot.data$stage,level=c("Egg","1-Cell","4-Cell","64-Cell","Sphere","Shield"))
filterfdr <- function(x){
     if (x < 0.05){return(-log10(x))}
     else{return(NA)}
}
filterpv <- function(x){
     if (x < 0.001){return(-log10(x))}
     else{return(NA)}
}
plot.data$fdr <- sapply(plot.data$fdr, filterfdr)
plot.data$p <- sapply(plot.data$p, filterpv)
#plot.data$fdr
#plot.data$term
#firebrick2
#darkblue
ggplot(plot.data, aes(x=stage, y=term)) + geom_point(aes(size = count,colour=p),na.rm=F)+scale_colour_gradient2(low='grey', high='darkblue', name='logp', na.value="grey10") + theme_bw() + xlab("") + ylab("")+
  theme(panel.background=element_blank(),
                panel.grid.major =element_blank(),
                panel.grid.minor = element_blank(),
                panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(family="Times",face="bold",size=10,color="black"),
        axis.text.y=element_text(family="Times",face="bold",size=10,color="black"),
        axis.title.x=element_text(family="Times",face="bold",size=18),
        axis.title.y=element_text(family="Times",face="bold",size=18),
        legend.title=element_text(family="Times",face="bold",size=16),
        legend.text=element_text(family="Times",size=16),
        strip.text.x=element_text(family="Times",size=16),
        strip.text.y=element_text(family="Times",size=16)) +
  scale_size_continuous(breaks=c(20,40,60),range=c(3,10))
  #scale_size_continuous(range=c(3,10))
#ggsave("top.fdr.pdf",width=8.3,height = 6)
ggsave("bottom.final.pdf", width=8, height = 5)
