script=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-new/code/bigMatrix.sh
matrixName="smooth-bigmatrix-stage6"
dir="/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10"
cd $dir
dir2=./
egg=egg-smooth.out
cell1=cell1-smooth.out
cell4=cell4-smooth.out
cell64=cell64-smooth.out
sphere=sphere-smooth.out
shield=shield-smooth.out
bash $script $dir2 $matrixName $egg $cell1 $cell4 $cell64 $sphere $shield
