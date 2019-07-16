
####################总的调用命令
cd /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23
##########1_step icshape replicate combine and make a huge matrix about icshape value
########## Notice to change the parameters in it like outpath
nohup bash ./code/command/step1.sh 


##2_step#################calculate diff value
###############all stage
script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/calculateDiff.R
cd w10_s10
mkdir 01_001_new
cd 01_001_new
nohup Rscript $script "../smooth-bigmatrix-stage6.txt"  0.01 ./ 


script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/calculateDiff.R
cd w10_s10
mkdir 01_005_new
cd 01_005_new
nohup Rscript $script "../smooth-bigmatrix-stage6.txt"  0.05 ./  &


cd /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10
mkdir 005_001_new
cp 01_001_new/res.txt 005_001_new/res.txt
mkdir 005_005_new
cp 01_005_new/res.txt 005_005_new/res.txt

script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/getDiff.R
cd /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/
Rscript $script res.txt  ./  0.1  0.01  0  & ##abs
Rscript $script res.txt  ./  0.1  0.01  -1 & ##down
Rscript $script res.txt  ./  0.1  0.01  1  &##up


cd /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_005_new/
Rscript $script res.txt  ./  0.1  0.05  0  & ##abs
Rscript $script res.txt  ./  0.1  0.05  -1 & ##down
Rscript $script res.txt  ./  0.1  0.05  1 & ##up


cd /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_001_new/
Rscript $script res.txt  ./  0.05  0.01  0  & ##abs
Rscript $script res.txt  ./  0.05  0.01  -1 & ##down
Rscript $script res.txt  ./  0.05  0.01  1 & ##up


cd /pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_005_new/
Rscript $script res.txt  ./  0.05  0.05  0   &##abs
Rscript $script res.txt  ./  0.05  0.05  -1 & ##down
Rscript $script res.txt  ./  0.05  0.05  1 & ##up



#3_step peak merge
script=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/code/command/step3.new.sh
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new
cd $dir
mkdir up down abs

for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/up  $dir/matrix.variable-up.txt $i
echo $i
done


for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/down $dir/matrix.variable-down.txt $i
echo $i
done


for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/abs  $dir/matrix.variable-abs.txt $i
echo $i
done


dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_005_new
cd $dir
mkdir up down abs

for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/up  $dir/matrix.variable-up.txt $i
echo $i
done


for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/down $dir/matrix.variable-down.txt $i
echo $i
done


for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/abs  $dir/matrix.variable-abs.txt $i
echo $i
done


dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_001_new
cd $dir
mkdir up down abs

for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/up  $dir/matrix.variable-up.txt $i
echo $i
done


for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/down $dir/matrix.variable-down.txt $i
echo $i
done


for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/abs  $dir/matrix.variable-abs.txt $i
echo $i
done



dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_005_new
cd $dir
mkdir up down abs

for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/up  $dir/matrix.variable-up.txt $i
echo $i
done


for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/down $dir/matrix.variable-down.txt $i
echo $i
done


for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
echo $i
bash $script  $dir/abs  $dir/matrix.variable-abs.txt $i
echo $i
done


#4_step Fimo找motif
####4_step FIMO
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new/
cd $dir
script=../../code/command/step4.new.sh
nohup bash $script $dir 

dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_005_new/
cd $dir 
script=../../code/command/step4.new.sh
nohup bash $script $dir 

dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_001_new/
cd  $dir 
script=../../code/command/step4.new.sh
nohup bash $script $dir 

dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/005_005_new/
cd  $dir 
script=../../code/command/step4.new.sh
nohup bash $script $dir 


#5_step Homer找富集的motif
############HOMER
cd $dir
script=../../code/command/step5.sh
nohup bash  $script $dir/abs 
nohup bash  $script $dir/up 
nohup bash  $script $dir/down 


###6_step 
#6_step ############################求相应的stable的集合
mkdir $dir/stable 
script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/getStable.R
Rscript $script $dir/stable ../../smooth-bigmatrix-stage6.txt  ../matrix.variable-up.txt ../matrix.variable-down.txt "egg_cell1:cell1_cell4:cell4_cell64:cell64_k1:k1_sphere:sphere_shield"

script=$dir/code/command/step3.sh
mkdir  stable/egg_cell1 stable/cell1_cell4 stable/cell4_cell64 stable/cell64_k1 stable/k1_sphere stable/sphere_shield

mv stable/egg_cell1-var.txt stable/egg_cell1
mv stable/cell1_cell4-var.txt stable/cell1_cell4
mv stable/cell4_cell64-var.txt stable/cell4_cell64
mv stable/cell64_k1-var.txt  stable/cell64_k1
mv stable/k1_sphere-var.txt  stable/k1_sphere
mv stable/sphere_shield-var.txt stable/sphere_shield

############peak merge
tranLen=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/GeneList/whole-trans/transcript-length.txt
intersectbed=/software/biosoft/software/MeRIP-PF/tools/BEDTools-Version-2.16.2/bin/intersectBed
path=$dir/stable

for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_k1 k1_sphere sphere_shield
do
	script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/seqMerge.pl 
	cd $path/$i
	perl $script 1 $i"-var.txt" $tranLen 20  20 50

	transBed=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/GeneList/whole-trans/wholeTransPositionDetail-40-200-40.bed
	$intersectbed -a window.used.bed -b $transBed  -f 0.6 -wa -wb  > window-anno.bed
	
	#######取序列做motif分析
	script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/seqMerge.pl 
	tranFa=/pnas/yangyg_group/zhangting/Reference/Zebrafish/Refseq_v10/transcriptome/transcriptome.fa
	bed=window-anno.bed
	perl $script 2 $bed $tranFa  > fg.fa
done














########################################下面的代码还没有调试


###############################各时期的stage window作为背景
#########调用/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/code/stage_all/01_001_new/Layout-new/Figure4A/getStable.R
########求各个时期stage的window

script=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/code/code/command/step3.sh
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/code/stage_all/01_001
cd $dir
mkdir  stable/egg_cell1 stable/cell1_cell4 stable/cell4_cell64 stable/cell64_k1 stable/k1_sphere stable/sphere_shield

mv stable/egg_cell1-var.txt stable/egg_cell1
mv stable/cell1_cell4-var.txt stable/cell1_cell4
mv stable/cell4_cell64-var.txt stable/cell4_cell64
mv stable/cell64_k1-var.txt  stable/cell64_k1
mv stable/k1_sphere-var.txt  stable/k1_sphere
mv stable/sphere_shield-var.txt stable/sphere_shield









































































































