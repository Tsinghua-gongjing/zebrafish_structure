
###inputfile=$1
###合并矩阵的名称=$2

# if [ $# -gt 0 ]; then
        # echo "参数个数为$#个"
# else
        # echo "没有参数"
# fi

cd ${1}


script=/pnas/yangyg_group/zhangting/software/icSHAPE-2/scripts/bigMatrix-correlated.pl
file=$3

i=1
for p in $*;
 do
	i=`expr $i + 1`
	if [ $i -gt 4 ]; then
	
		file=$file:$p
		
	fi
 done;
 
echo $file
perl $script -i $file -o ${2}.txt

















