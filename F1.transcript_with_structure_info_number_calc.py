import os
import sys
from glob import glob
def readicshapeout(fname, p):
    tx_dict = {}
    with open(fname, 'r') as f:
        for l in f:
            splitLine = l.strip().split("\t")
            tx_id = splitLine[0]
            length = int(splitLine[1])
            scores = splitLine[3:]
            assert length == len(scores)
            if length - scores.count('NULL') > 0:
                tx_dict[tx_id] = 1
    return tx_dict


def main():
    wd = '/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win'
    f_ls = glob(wd+"/*icshape.w200.s30.T2.t200.out")
    for f in f_ls:
        print "\t".join([f,str(len(readicshapeout(f, float(sys.argv[1]))))])    


if __name__ == "__main__":
    main()



# /Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/egg.icshape.w200.s30.T2.t200.out 6070
# /Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/1cell.icshape.w200.s30.T2.t200.out       6004
# /Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/4cell.icshape.w200.s30.T2.t200.out       6910
# /Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/64cell.icshape.w200.s30.T2.t200.out      7491
# /Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out      7655
# /Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out      6187