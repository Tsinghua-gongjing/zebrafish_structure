export PATH=$PATH:/pnas/yangyg_group/yangxin/software/Homer/bin
shufflebed=/software/biosoft/software/MeRIP-PF/tools/BEDTools-Version-2.16.2/bin/shuffleBed
tranSize=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/GeneList/whole-trans/transcript-length.txt
tranfa=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/transcriptome/transcriptome.fa
known_motif=/pnas/yangyg_group/zhangting/SBY_strucutre-m5C-layout/supplementary/moleculalr_cell/HOMER/motif/rna.motifs

dir=$1
for i in "egg_cell1"  "cell1_cell4" "cell4_cell64" "cell64_k1" "k1_sphere" "sphere_shield"
do
	cd $dir/$i
	# $shufflebed -i ./window-anno2.bed -chrom -g $tranSize > shuffle_window-anno.bed
	# $shufflebed -i ./window-anno2_utr5.bed -chrom -g $tranSize > shuffle_window-anno_utr5.bed
	# $shufflebed -i ./window-anno2_cds.bed -chrom -g $tranSize  > shuffle_window-anno_cds.bed
	# $shufflebed -i ./window-anno2_utr3.bed -chrom -g $tranSize  > shuffle_window-anno_utr3.bed
	
	$shufflebed -i ./window-anno.bed -g $tranSize > shuffle_window-anno2.bed
	$shufflebed -i ./window-anno_utr5.bed -g $tranSize > shuffle_window-anno2_utr5.bed
	$shufflebed -i ./window-anno_cds.bed -g $tranSize  > shuffle_window-anno2_cds.bed
	$shufflebed -i ./window-anno_utr3.bed -g $tranSize  > shuffle_window-anno2_utr3.bed
	
	
	
	mkdir -p Homer2/all Homer2/utr5 Homer2/cds Homer2/utr3
	
	######################all
	findMotifsGenome.pl ./window-anno.bed $tranfa ./Homer2/all/Annotation_5nt/ -size given -p 1 -len 5 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno.bed $tranfa ./Homer2/all/Annotation_6nt/ -size given -p 1 -len 6 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno.bed $tranfa ./Homer2/all/Annotation_7nt/ -size given -p 1 -len 7 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno.bed $tranfa ./Homer2/all/Annotation_8nt/ -size given -p 1 -len 8 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno.bed $tranfa ./Homer2/all/Annotation_9nt/ -size given -p 1 -len 9 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2.bed -mcheck $known_motif -mknown $known_motif &

	#####################utr5
	findMotifsGenome.pl ./window-anno_utr5.bed $tranfa ./Homer2/utr5/Annotation_5nt/ -size given -p 1 -len 5 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_utr5.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_utr5.bed $tranfa ./Homer2/utr5/Annotation_6nt/ -size given -p 1 -len 6 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_utr5.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_utr5.bed $tranfa ./Homer2/utr5/Annotation_7nt/ -size given -p 1 -len 7 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_utr5.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_utr5.bed $tranfa ./Homer2/utr5/Annotation_8nt/ -size given -p 1 -len 8 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_utr5.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_utr5.bed $tranfa ./Homer2/utr5/Annotation_9nt/ -size given -p 1 -len 9 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_utr5.bed -mcheck $known_motif -mknown $known_motif &
	

	#####################cds
	findMotifsGenome.pl ./window-anno_cds.bed $tranfa ./Homer2/cds/Annotation_5nt/ -size given -p 1 -len 5 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_cds.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_cds.bed $tranfa ./Homer2/cds/Annotation_6nt/ -size given -p 1 -len 6 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_cds.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_cds.bed $tranfa ./Homer2/cds/Annotation_7nt/ -size given -p 1 -len 7 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_cds.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_cds.bed $tranfa ./Homer2/cds/Annotation_8nt/ -size given -p 1 -len 8 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_cds.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_cds.bed $tranfa ./Homer2/cds/Annotation_9nt/ -size given -p 1 -len 9 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_cds.bed -mcheck $known_motif -mknown $known_motif &
	
	
	############utr3
	findMotifsGenome.pl ./window-anno_utr3.bed $tranfa ./Homer2/utr3/Annotation_5nt/ -size given -p 1 -len 5 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_utr3.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_utr3.bed $tranfa ./Homer2/utr3/Annotation_6nt/ -size given -p 1 -len 6 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_utr3.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_utr3.bed $tranfa ./Homer2/utr3/Annotation_7nt/ -size given -p 1 -len 7 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_utr3.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_utr3.bed $tranfa ./Homer2/utr3/Annotation_8nt/ -size given -p 1 -len 8 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_utr3.bed -mcheck $known_motif -mknown $known_motif &

	findMotifsGenome.pl ./window-anno_utr3.bed $tranfa ./Homer2/utr3/Annotation_9nt/ -size given -p 1 -len 9 -rna -chopify -norevopp -cache 1000 -bg shuffle_window-anno2_utr3.bed -mcheck $known_motif -mknown $known_motif &
	
done


 



