#rna2meme 
function FIMO(){
	cd $1
	fimo=~/bin/fimo
	fg=fg.fa
	fg_utr5=fg_utr5.fa
	fg_cds=fg_cds.fa
	fg_utr3=fg_utr3.fa

	motif=/pnas/yangyg_group/zhangting/software/meme_4.12.0/motif/motif_CISBP_RNA_narrow/Collapsed.meme
	
	nohup $fimo --oc fimo_out-all  --thresh 0.001 $motif  $fg &
	nohup $fimo --oc fimo_out-utr5  --thresh 0.001 $motif  $fg_utr5 &
	nohup $fimo --oc fimo_out-cds  --thresh 0.001 $motif  $fg_cds &
	nohup $fimo --oc fimo_out-utr3  --thresh 0.001 $motif  $fg_utr3 &

}

for j in "egg_cell1"  "cell1_cell4" "cell4_cell64" "cell64_sphere" "sphere_shield"
do
		dir=$1/up/$j
		echo $dir
		FIMO  $dir &
done


for j in "egg_cell1"  "cell1_cell4" "cell4_cell64" "cell64_sphere" "sphere_shield"
do
		dir=$1/down/$j
		echo $dir
		FIMO  $dir &
done


for j in "egg_cell1"  "cell1_cell4" "cell4_cell64" "cell64_sphere" "sphere_shield"
do
		dir=$1/abs/$j
		echo $dir
		FIMO  $dir &
done


#for j in "egg_cell1"  "cell1_cell4" "cell4_cell64" "cell64_sphere" "sphere_shield"
#do
#		dir=$1/stable/$j
#		echo $dir
#		FIMO  $dir &
		
#done

	
	


	




	
	
	
	
	
	
