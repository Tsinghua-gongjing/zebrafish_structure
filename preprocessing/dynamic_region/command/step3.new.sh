##############差异window上下左右外扩20base，然后将相邻的window合并打断成50nt的peak
####
####$1=dir
###$2=inputFile
	script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/seqMerge.pl 
	tranLen=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/GeneList/whole-trans/transcript-length.txt
	intersectbed=/software/biosoft/software/MeRIP-PF/tools/BEDTools-Version-2.16.2/bin/intersectBed
	echo $1
	echo $2
	cd $1
	
	matrix=$2
	/bin/mkdir egg_cell1 cell1_cell4 cell4_cell64 cell64_sphere sphere_shield
	awk ' $2!="NA"' $matrix >  egg_cell1/egg_cell1-var.txt
	awk ' $3!="NA"' $matrix >   cell1_cell4/cell1_cell4-var.txt
	awk ' $4!="NA" '  $matrix > cell4_cell64/cell4_cell64-var.txt
	awk ' $5!="NA" '  $matrix > cell64_sphere/cell64_sphere-var.txt
	awk  ' $6!="NA" '  $matrix > sphere_shield/sphere_shield-var.txt
	
	outputPath=$3
	cd $outputPath
	
	perl $script 1 $outputPath"-var.txt" $tranLen 10  10 30

	transBed=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/GeneList/whole-trans/wholeTransPositionDetail.bed
	$intersectbed -a window.used.bed -b $transBed -wa -wb  > window-anno.bed

	#########区分utr5 cds utr3
	grep "utr5" window-anno.bed> window-anno_utr5.bed
	grep "cds" window-anno.bed> window-anno_cds.bed
	grep "utr3" window-anno.bed> window-anno_utr3.bed
	
	#######取序列做motif分析
	script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/seqMerge.pl 
	tranFa=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/transcriptome/transcriptome.fa
	
	perl $script 2 window-anno.bed $tranFa  > fg.fa
	perl $script 2 window-anno_utr5.bed $tranFa  > fg_utr5.fa
	perl $script 2 window-anno_cds.bed $tranFa  > fg_cds.fa
	perl $script 2 window-anno_utr3.bed $tranFa  > fg_utr3.fa


	# maternal=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_maternal/maternal-gene/DMSO/maternal-decay.txt
	# maternal_Decay=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_maternal/maternal-gene/FC1.5/MD.txt
	# zygotic_Decay=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_maternal/maternal-gene/FC1.5/ZD.txt
	
	# awk -v OFS="\t" 'NR==FNR{x[$1]=2}NR>FNR{if(x[$1]){print $0}}' $maternal window-anno.bed > maternal_window.anno.bed
	# awk -v OFS="\t" 'NR==FNR{x[$1]=2}NR>FNR{if(x[$1]){print $0}}' $maternal_Decay window-anno.bed > maternal_Decay_anno.used.bed
	# awk -v OFS="\t" 'NR==FNR{x[$1]=2}NR>FNR{if(x[$1]){print $0}}' $zygotic_Decay window-anno.bed > zygotic_Decay_anno.used.bed


	
	#script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/seqMerge.pl 
	#tranFa=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/transcriptome/transcriptome.fa
	
	#bed=maternal_window.anno.bed
	#perl $script 2 $bed $tranFa  > maternal-fg.fa

	# bed=maternal_Decay_anno.used.bed
	# perl $script 2 $bed $tranFa  > maternal-Decay-fg.fa

	# bed=zygotic_Decay_anno.used.bed
	# perl $script 2 $bed $tranFa  >  zygotic-Decay-fg.fa





