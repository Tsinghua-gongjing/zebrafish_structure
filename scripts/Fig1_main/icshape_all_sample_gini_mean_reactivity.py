import subprocess
import itertools
from nested_dict import nested_dict
import sys,os
import pandas as pd, numpy as np
from multiprocessing import Pool
import gj

def readIc(ic_file):
    "Read icSHAPE Reactivity File"
    icDict = dict()
    IN = open(ic_file)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        icDict[arr[0]] = arr[3:]
        line = IN.readline()
    IN.close()
    return icDict

def loadTransGtfBed2(ref_bed='/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/ref_GRCz10_top_level.trans.bed.2'):
    H = open(ref_bed)
    line = H.readline()
    trans_dict = nested_dict()
    header_ls = ['tx', 'gene', 'type', 'length', 'utr_5_start', 'utr_5_end', 'cds_start', 'cds_end', 'utr_3_start', 'utr_3_end']
    while line:
        if line.startswith('#'): line = H.readline(); continue
        arr = line.strip().split()
        for i,j in zip(header_ls, arr):
            trans_dict[arr[0]][i] = j
        line = H.readline()
    H.close()
    print "read: %s, n=%s"%(ref_bed, len(trans_dict))
    return trans_dict.to_dict()

def gini(list_of_values,mode='mean_reactivity',null_pct=1):
    if len(list_of_values) == 0: return -1
    if list_of_values.count('NULL')/float(len(list_of_values)) > null_pct: return -1
    list_of_values = [i for i in list_of_values if i != 'NULL']
    if len(list_of_values) == 0: return -1
    if type(list_of_values[0]) is str:
        list_of_values = map(float,list_of_values)
    if mode == 'mean_reactivity':
        return np.mean(map(float,list_of_values))
    if mode == 'gini':
        if sum(list_of_values) == 0: return 0.67
        sorted_list = sorted(list_of_values)
        height, area = 0, 0
        for value in sorted_list:
            height += value
            area += height - value / 2.
        fair_area = height * len(list_of_values) / 2.
        return (fair_area - area) / fair_area

def all_gini(RPKM_combine=None, mode='gini', null_pct=1):
    if RPKM_combine is None:
        RPKM_combine = '/Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelationPairwiseNew/2018_03_12_gini/RPKM_combine.merge.txt'
    df = pd.read_csv(RPKM_combine, header=None, index_col=0, sep='\t')
    print df.head()

    trans_dict = loadTransGtfBed2('/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.trans.bed2')

    sample_ls = ['egg', '1cell', '4cell', '64cell', '1K', 'sphere', 'shield']
    sample_path = ['/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/%s.icshape.w200.s30.T2.t200.out'%(i) for i in sample_ls] # norm by window

    sample_ic_dict =  nested_dict()
    for i,j in zip(sample_path, sample_ls):
        print "read icshape: %s"%(i)
        sample_ic_dict[j] = readIc(i)

    t_cutoff = sample_path[0].split('.')[-2]
    savefn = RPKM_combine.replace('.txt', '.%s.%s.null%s.txt'%(t_cutoff, mode, int(null_pct*100)))
    
    SAVEFN = open(savefn, 'w')
    print >>SAVEFN, '\t'.join(['%s(transcript)\t%s(UTR5)\t%s(CDS)\t%s(UTR3)'%(i,i,i,i) for i in sample_ls])
    for tx in df.index:
        gini_ls = [tx]
        for i in sample_ls:
            if sample_ic_dict[i].has_key(tx) and trans_dict.has_key(tx):
                utr_5_start, utr_5_end, cds_start, cds_end, utr_3_start, utr_3_end = [int(trans_dict[tx][g]) for g in ['utr_5_start', 'utr_5_end', 'cds_start', 'cds_end', 'utr_3_start', 'utr_3_end']] 
                if utr_5_start == 0:
                    utr_5_gini = 'NULL'
                else:
                    utr_5_gini = gini(sample_ic_dict[i][tx][(utr_5_start-1):(utr_5_end)], mode=mode, null_pct=null_pct)
                    if utr_5_gini < 0:
                        utr_5_gini = 'NULL'
                if utr_3_start == 0:
                    utr_3_gini = 'NULL'
                else:
                    utr_3_gini = gini(sample_ic_dict[i][tx][(utr_3_start-1):(utr_3_end)], mode=mode, null_pct=null_pct)
                    if utr_3_gini < 0:
                        utr_3_gini = 'NULL'
                cds_gini = gini(sample_ic_dict[i][tx][(cds_start-1):(cds_end)], mode=mode, null_pct=null_pct)
                if cds_gini < 0:
                    cds_gini = 'NULL'
                transcript_gini = gini(sample_ic_dict[i][tx][0:], mode=mode, null_pct=null_pct)
                if transcript_gini < 0:
                    transcript_gini = 'NULL'
                sample_gini_ls = [transcript_gini, utr_5_gini, cds_gini, utr_3_gini]
            else:
                sample_gini_ls = ['NULL','NULL', 'NULL', 'NULL']
            gini_ls += sample_gini_ls
        print >> SAVEFN, '\t'.join(map(str, gini_ls))
    SAVEFN.close()

def all_gini_run():
    for mode in ['gini', 'mean_reactivity']:
        for null_pct in [0.4]:
            all_gini(RPKM_combine=None, mode=mode, null_pct=null_pct)

def main():
    all_gini_run()

if __name__ == '__main__':
    main()