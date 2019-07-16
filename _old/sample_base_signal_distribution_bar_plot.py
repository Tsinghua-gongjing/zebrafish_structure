import gj
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
import sys
from nested_dict import nested_dict
import pandas as pd
import numpy as np
from pyfasta import Fasta
import structure
from multiprocessing import Pool
import matplotlib.backends.backend_pdf
import scipy.stats as stats

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

def read_compare(txt=None):
    distance_ls = []
    base_ls = []
    with open(txt, 'r') as TXT:
        for n, line in enumerate(TXT):
            if n % 1000000 == 0:
                print "read line: %s" % (n)
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            arr = line.split('\t')
            if arr[3] == 'NULL' or arr[4] == 'NULL':
                continue
            else:
                distance = float(arr[4]) - float(arr[3])
                distance_ls.append(distance)
                base_ls.append(arr[2].upper())
    return distance_ls, base_ls

def read_fa(fa='/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.fa'):
    fa_dict1 = Fasta(fa, key_fn=lambda key: key.split("\t")[0])
    fa_dict = {i.split()[0]: j[0:] for i, j in fa_dict1.items()}
    print fa_dict.keys()[0:3]
    return fa_dict

def opening_base_compare_stats_bar_plot():
    egg_cell1 = '/Share/home/zhangqf7/gongjing/zebrafish/result/structure_change_compare/egg_1cell.base.txt'
    sphere_shield = '/Share/home/zhangqf7/gongjing/zebrafish/result/structure_change_compare/sphere_shield.base.txt'
    egg_cell1_distance_ls, egg_cell1_base_ls = read_compare(egg_cell1)
    sphere_shield_distance_ls, sphere_shield_base_ls = read_compare(
        sphere_shield)
    print gj.printFuncRun('creat df_egg_cell1')
    df_egg_cell1 = pd.DataFrame({'difference': egg_cell1_distance_ls,
                       'sample': '1cell-egg', 'base': egg_cell1_base_ls})
    print gj.printFuncRun('finsh creat df_egg_cell1')
    print gj.printFuncRun('creat df_sphere_shield')
    df_sphere_shield = pd.DataFrame({'difference': sphere_shield_distance_ls,
                       'sample': 'shield-sphere', 'base': sphere_shield_base_ls})
    print gj.printFuncRun('finish creat df_sphere_shield')
    # print gj.printFuncRun('concat')
    # df = pd.concat([df_egg_cell1, df_sphere_shield], axis=0)
    # print gj.printFuncRun('concat')
    # print df.head()

    diff_cutoff = 0.25
    df_egg_cell1_above = df_egg_cell1[df_egg_cell1['difference']>=diff_cutoff]
    df_sphere_shield_above = df_sphere_shield[df_sphere_shield['difference']>=diff_cutoff]
    print 'egg-1cell', 'all', df_egg_cell1['base'].value_counts(), df_egg_cell1_above['base'].value_counts()
    print "sphere-shield", 'all', df_sphere_shield['base'].value_counts(), df_sphere_shield_above['base'].value_counts()
    egg_1cell_all = df_egg_cell1['base'].value_counts().to_dict()
    egg_1cell_above = df_egg_cell1_above['base'].value_counts().to_dict()
    sphere_shield_all = df_sphere_shield['base'].value_counts().to_dict()
    sphere_shield_above = df_sphere_shield_above['base'].value_counts().to_dict()

    base_ls, sample_ls, value_ls = [],[],[]
    # for d,label in zip([egg_1cell_all, egg_1cell_above, sphere_shield_all, sphere_shield_above], ['egg_1cell_all', 'egg_1cell_above', 'sphere_shield_all', 'sphere_shield_above']):
    for d,label in zip([egg_1cell_all, egg_1cell_above], ['egg_1cell_all', 'egg_1cell_above']):
        if 'N' in d: d.pop('N')
        for b in ['A', 'T', 'C', 'G']:
            base_ls.append(b)
            sample_ls.append(label)
            value = d[b]
            value_ls.append(value)

            value_ratio = d[b] / float(sum(d.values()))
            base_ls.append(b)
            sample_ls.append(label+'\n(ratio)')
            value_ls.append(value_ratio)
    df_stat = pd.DataFrame({'base':base_ls, 'sample':sample_ls, 'value':value_ls})
    print df_stat

    fig,ax = plt.subplots(2,1)
    sns.barplot(x='sample', y='value', hue='base', data=df_stat[df_stat['value']>=1], hue_order=['A', 'T', 'C', 'G'], ax=ax[0])
    sns.barplot(x='sample', y='value', hue='base', data=df_stat[df_stat['value']<1], hue_order=['A', 'T', 'C', 'G'], ax=ax[1])
    plt.tight_layout()
    savefn = '/Share/home/zhangqf7/gongjing/zebrafish/result/structure_change_compare/base_num_ratio.onlyegg_1cell.pdf'
    plt.savefig(savefn)
    plt.close()

    print "all",egg_1cell_all,"above",egg_1cell_above
    for b in ['A', 'T', 'C', 'G']:
        b_above = egg_1cell_above[b]
        non_b_above = sum(egg_1cell_above.values()) - b_above
        b_stable = egg_1cell_all[b] - b_above
        non_b_stable = sum(egg_1cell_all.values()) - b_above - non_b_above - b_stable
        print b_above, non_b_above, b_stable, non_b_stable
        oddsratio, pvalue = stats.fisher_exact([[b_above, non_b_above], [b_stable, non_b_stable]])
        print b, oddsratio, pvalue

def main():
    opening_base_compare_stats_bar_plot()

if __name__ == '__main__':
    main()