path1=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_RBP/window/t200/abs/01/
path2=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_RBP/window/t200/abs/01_005/
path3=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_RBP/window/t200/down/01/
path4=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_RBP/window/t200/down/01_005/
path5=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_RBP/window/t200/up/01/
path6=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_RBP/window/t200/up/01_005/

script=/pnas/yangyg_group/zhangting/SBY_structure-m5C-new/0-analysis/0_RBP/code/RBP-staticstics.R
nohup Rscript  $script  $path1 "utr3-motif-d/fimo_out-all"  "abs-01-all"  0.01 &
nohup Rscript  $script  $path1 "utr3-motif-d/fimo_out-maternal"  "abs-01-maternal" 0.01 &
nohup Rscript  $script  $path1 "utr3-motif-d/fimo_out-md"  "abs-01-md" 0.01 &
nohup Rscript  $script  $path1 "utr3-motif-d/fimo_out-zd"  "abs-01-zd" 0.01 &


nohup Rscript  $script  $path2 "utr3-motif-d/fimo_out-all"  "abs-01-005-all"  0.01 &
nohup Rscript  $script  $path2 "utr3-motif-d/fimo_out-maternal"  "abs-01-005-maternal" 0.01 &
nohup Rscript  $script  $path2 "utr3-motif-d/fimo_out-md"  "abs-01-005-md" 0.01 &
nohup Rscript  $script  $path2 "utr3-motif-d/fimo_out-zd"  "abs-01-005-zd" 0.01 &

nohup Rscript  $script  $path3 "utr3-motif-d/fimo_out-all"  "down-01-all"  0.01 &
nohup Rscript  $script  $path3 "utr3-motif-d/fimo_out-maternal"  "down-01-maternal" 0.01 &
nohup Rscript  $script  $path3 "utr3-motif-d/fimo_out-md"  "down-01-md" 0.01 &
nohup Rscript  $script  $path3 "utr3-motif-d/fimo_out-zd"  "down-01-zd" 0.01 &

nohup Rscript  $script  $path4 "utr3-motif-d/fimo_out-all"  "down-01-005-all"  0.01 &
nohup Rscript  $script  $path4 "utr3-motif-d/fimo_out-maternal"  "down-01-005-maternal" 0.01 &
nohup Rscript  $script  $path4 "utr3-motif-d/fimo_out-md"  "down-01-005-md" 0.01 &
nohup Rscript  $script  $path4 "utr3-motif-d/fimo_out-zd"  "down-01-005-zd" 0.01 &

nohup Rscript  $script  $path5 "utr3-motif-d/fimo_out-all"  "up-01-all"  0.01 &
nohup Rscript  $script  $path5 "utr3-motif-d/fimo_out-maternal"  "up-01-maternal" 0.01 &
nohup Rscript  $script  $path5 "utr3-motif-d/fimo_out-md"  "up-01-md" 0.01 &
nohup Rscript  $script  $path5 "utr3-motif-d/fimo_out-zd"  "up-01-zd" 0.01 &


nohup Rscript  $script  $path6 "utr3-motif-d/fimo_out-all"  "up-01-005-all"  0.01 &
nohup Rscript  $script  $path6 "utr3-motif-d/fimo_out-maternal"  "up-01-005-maternal" 0.01 &
nohup Rscript  $script  $path6 "utr3-motif-d/fimo_out-md"  "up-01-005-md" 0.01 &
nohup Rscript  $script  $path6 "utr3-motif-d/fimo_out-zd"  "up-01-005-zd" 0.01 &









 







