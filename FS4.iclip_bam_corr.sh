# https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html#example

export PYTHONPATH=/Share/home/zhangqf7/usr/numpy-1.16.0/lib/python2.7/site-packages/

gene_bed=/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/danRerV10.refGene.genes.bed

# used sam need sorted & indexed, otherwise error
rep1=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-4/bwa/iCLIP.sorted.bam
rep2=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-10/bwa/iCLIP.sorted.bam
rep3=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-16/bwa/iCLIP.sorted.bam
rep4=/Share/home/zhangqf7/jinsong_zhang/zebrafish/data/iclip/20181224/Rawdata/shi-zi-22/bwa/iCLIP.sorted.bam

sample=h6_low
out_npz=/Share/home/zhangqf7/gongjing/zebrafish/result/iCLIP/bam_corr/${sample}.readCounts.gene.npz
out_tab=/Share/home/zhangqf7/gongjing/zebrafish/result/iCLIP/bam_corr/${sample}.readCounts.gene.tab

multiBamSummary bins \
		  --bamfiles $rep1 $rep2 $rep3 $rep4 \
		  --labels rep1 rep2 rep3 rep4 \
		  -out $out_npz \
		  --outRawCounts $out_tab \
		  --BED $gene_bed

out_png=/Share/home/zhangqf7/gongjing/zebrafish/result/iCLIP/bam_corr/${sample}.corr.gene.png

plotCorrelation \
		  -in $out_npz \
		  --corMethod spearman --skipZeros \
		  --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
		  -o $out_png