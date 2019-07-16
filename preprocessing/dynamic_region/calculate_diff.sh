script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/calculateDiff.R
outdir="/pnas/yangyg_group/zhangting/SBY_structure-m5C-new2/0_RBP/window-may23/w10_s10"
inputfile="smooth-bigmatrix-stage6.txt"
cd $outdir
mkdir 01_001_new
cd 01_001_new
Rscript $script $outdir/$inputfile 0.01 ./

script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/getDiff.R
cd $outdir/01_001_new
Rscript $script res.txt  ./  0.1  0.01  0  ##abs
Rscript $script res.txt  ./  0.1  0.01  -1 ##down
Rscript $script res.txt  ./  0.1  0.01  1 ##up

