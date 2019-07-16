script=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/code/command/step3.new.sh
dir=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10/01_001_new
cd $dir
mkdir up down abs
for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_sphere sphere_shield
do
echo $i
bash $script $dir/up $dir/matrix.variable-up.txt $i
echo $i
done

for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_sphere sphere_shield
do
echo $i
bash $script $dir/down $dir/matrix.variable-down.txt $i
echo $i
done

for i in egg_cell1 cell1_cell4 cell4_cell64 cell64_sphere sphere_shield
do
echo $i
bash $script $dir/abs $dir/matrix.variable-abs.txt $i
echo $i
done

