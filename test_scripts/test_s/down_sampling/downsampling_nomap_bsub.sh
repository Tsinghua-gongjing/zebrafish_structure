
for n in {1..100}
do

bsub -q Z-ZQF -n 2 -eo /Share/home/zhangqf7/gongjing/zebrafish/script/downsampling/$1_down${n}.err -oo /Share/home/zhangqf7/gongjing/zebrafish/script/downsampling/$1_down${n}.out bash downsampling_nomap.sh $1 $2 ${n}

done
