#!/bin/sh
#PBS -q core24
#PBS -l mem=50gb,walltime=500:00:00,nodes=1:ppn=14
#PBS -o icshape_stdout-1.txt
#PBS -e icshape_stderr-1.txt
#HSCHED -s icshape+tophate+zebrafishV10

<<!
#####################将两个批次的数据combine在一起
transcript=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/transcriptome/transcriptome
script=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_RBP/code/icshape-combine/icshape_combine.sh

######egg
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_combine/egg/
dmso1=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-8/rt/shi-zf-8.combined.transcriptome.rt
dmso2=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-24/rt/shi-zf-24.combined.transcriptome.rt
nai1=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-9/rt/shi-zf-9.combined.transcriptome.rt
nai2=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-25/rt/shi-zf-25.combined.transcriptome.rt
nohup bash $script $transcript $dir $dmso1 $dmso2 $nai1 $nai2 &


######cell1
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_combine/cell1/
dmso3=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-10/rt/shi-zf-10.combined.transcriptome.rt
dmso4=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-26/rt/shi-zf-26.combined.transcriptome.rt
nai3=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-11/rt/shi-zf-11.combined.transcriptome.rt
nai4=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-27/rt/shi-zf-27.combined.transcriptome.rt
nohup bash $script $transcript $dir $dmso3 $dmso4 $nai3 $nai4 &


######cell4
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_combine/cell4/
dmso5=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-12/rt/shi-zf-12.combined.transcriptome.rt
dmso6=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-28/rt/shi-zf-28.combined.transcriptome.rt
nai5=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-13/rt/shi-zf-13.combined.transcriptome.rt
nai6=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-29/rt/shi-zf-29.combined.transcriptome.rt
nohup  bash $script $transcript $dir $dmso5 $dmso6 $nai5 $nai6 &

######cell64
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_combine/cell64/
dmso7=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-14/rt/shi-zf-14.combined.transcriptome.rt
dmso8=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-30/rt/shi-zf-30.combined.transcriptome.rt
nai7=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-15/rt/shi-zf-15.combined.transcriptome.rt
nai8=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-31/rt/shi-zf-31.combined.transcriptome.rt
nohup  bash $script $transcript $dir $dmso7 $dmso8 $nai7 $nai8 &

#####1k
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_combine/k1/
dmso9=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-16/rt/shi-zf-16.combined.transcriptome.rt
dmso10=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-32/rt/shi-zf-32.combined.transcriptome.rt
nai9=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-20/rt/shi-zf-20.combined.transcriptome.rt
nai10=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-33/rt/shi-zf-33.combined.transcriptome.rt
nohup  bash $script $transcript $dir $dmso9 $dmso10 $nai9 $nai10 &

###sphere
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_combine/sphere/
dmso11=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-17/rt/shi-zf-17.combined.transcriptome.rt
dmso12=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-40/rt/shi-zf-40.combined.transcriptome.rt
nai11=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-21/rt/shi-zf-21.combined.transcriptome.rt
nai12=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-35/rt/shi-zf-35.combined.transcriptome.rt
nohup  bash $script $transcript  $dir $dmso11 $dmso12 $nai11 $nai12 &

###shield
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_combine/shield/
dmso13=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-18/rt/shi-zf-18.combined.transcriptome.rt
dmso14=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-36/rt/shi-zf-36.combined.transcriptome.rt
nai13=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-22/rt/shi-zf-22.combined.transcriptome.rt
nai14=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-37/rt/shi-zf-37.combined.transcriptome.rt
nohup  bash $script $transcript  $dir $dmso13 $dmso14 $nai13 $nai14 &


###epiboly
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_combine/epiboly/
dmso15=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-19/rt/shi-zf-19.combined.transcriptome.rt
dmso16=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-38/rt/shi-zf-38.combined.transcriptome.rt
nai15=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-23/rt/shi-zf-23.combined.transcriptome.rt
nai16=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/icshape-zf-39/rt/shi-zf-39.combined.transcriptome.rt
nohup  bash $script $transcript  $dir $dmso15 $dmso16 $nai15 $nai16 &
!

function Calculation(){
	script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/icshape-smooth.sh
	positionFile=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/GeneList/whole-trans/wholeTransPositionDetail.txt
	#####################smooth
	window=$1
	step=$2
	outPath=$3
	##########create smoothed window
	cd $outPath
	mkdir "rep_w"$1"_s"$2
	dir=$outPath/"rep_w"$1"_s"$2
	
	
	input1=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_combine/rep1/egg.icshape.w200.s30.T2.t200.out
	output='rep1-egg-smooth.out'
	bash $script $dir $output $window $step $input1 $positionFile

	input2=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_combine/rep2/egg.icshape.w200.s30.T2.t200.out
	output='rep2-egg-smooth.out'
	bash $script $dir $output $window $step $input2 $positionFile

	input3=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_combine/rep1/1cell.icshape.w200.s30.T2.t200.out
	output='rep1-cell1-smooth.out'
	bash $script $dir $output $window $step $input1 $positionFile

	input4=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_combine/rep2/1cell.icshape.w200.s30.T2.t200.out
	output='rep2-cell1-smooth.out'
	bash $script $dir $output $window $step $input2 $positionFile

	input5=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_combine/rep1/4cell.icshape.w200.s30.T2.t200.out
	output='rep1-cell4-smooth.out'
	bash $script $dir $output $window $step $input1 $positionFile

	input6=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_combine/rep2/4cell.icshape.w200.s30.T2.t200.out
	output='rep2-cell4-smooth.out'
	bash $script $dir $output $window $step $input2 $positionFile
	
	#####################合成大矩阵并求差异的window
	###outputDir=$1
	###合并矩阵的名称=$2
	###需要合并的文件:$3....
	script=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-new/code/bigMatrix.sh
	matrixName="smooth-bigmatrix-stage6"

	###########$dirName
	cd $dir
	dir2=./
	input1=rep1-egg-smooth.out
	input2=rep2-egg-smooth.out
	input3=rep1-cell1-smooth.out
	input4=rep2-cell1-smooth.out
	input5=rep1-cell4-smooth.out
	input6=rep2-cell4-smooth.out
	bash $script $dir2 $matrixName $input1  $input2 $input3 $input4 $input5 $input6
}

Calculation 10 10 /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23


