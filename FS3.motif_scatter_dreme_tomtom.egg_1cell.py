import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Helvetica"
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
from nested_dict import nested_dict
import sys,os
import pandas as pd, numpy as np
from scipy import stats
from adjustText import adjust_text


def read_RBP_name():
    txt = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/dynamic_region/RBP.txt'
    df = pd.read_csv(txt, header=None, sep='\t')
    return {i:j for i,j in zip(df[1], df[2])}

def read_dmeme(dreme):
    import re
    motif_ls = []
    with open(dreme, 'r') as TXT:
        for line in TXT:
            line = line.strip()
            arr = re.split('\s+', line)
            if 'BEST' not in line:
                continue
            print arr
            motif_ls.append(arr)
    motif_df = pd.DataFrame(motif_ls)
    return motif_df

def read_tomtom(tomtom, rbp_name):
    df = pd.read_csv(tomtom, header=0, sep='\t')
    print df.head()
    d = nested_dict(1, list)
    for motif,target,pvalue in zip(df['Query_ID'], df['Target_ID'], df['p-value']):
        if rbp_name.has_key(target):
            target = rbp_name[target]
        d[motif].append('%s(%s)'%(target, pvalue))
    return d.to_dict()

def plot_motif_logo(motif_ls, site_num_ls, pvalue_ls, tomtom_file, rbp_name):
    dd = pd.DataFrame({'site_num':site_num_ls, 'pvalue':pvalue_ls})
    dd['-log10(pvalue)'] = -np.log10(dd['pvalue'])
    dd['log2(site_num)'] = np.log2(dd['site_num'])
    dd['motif'] = motif_ls
    print dd

    sns.set_style('ticks')
    fig,ax = plt.subplots(figsize=(25,20))
    dd.plot(kind='scatter', x='log2(site_num)', y='-log10(pvalue)', ax=ax)
    ax.set_xlim(0,)
    
    tomtom = read_tomtom(tomtom_file, rbp_name)
    
    texts = []
    for x,y,t in zip(dd['log2(site_num)'], dd['-log10(pvalue)'], dd['motif']):
        if tomtom.has_key(t):
            RBP = t + '\n'+'\n'.join(tomtom[t])
        else:
            RBP = t
        if t in ['TDTWV', 'TDWDT', 'AMACAMAC', 'TAANC', 'ACAYT']:
            texts.append(plt.text(x, y, RBP, fontsize=12))
    adjust_text(texts, only_move={'text': 'xy'})
    
    plt.title(tomtom_file.split('/')[-2])
    
    plt.tight_layout()
    savefn = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/FS3.motif_scatter_dreme_tomtom.egg_1cell.pdf'
    plt.savefig(savefn)

rbp_name = read_RBP_name()
print rbp_name

motif_df = read_dmeme('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/dynamic_region/egg_cell1/utr3/dreme_egg_1cell.txt')
motif_df.to_csv('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/dynamic_region/egg_cell1/utr3/egg_cell1_motif_search_site_vs_pval.txt', header=True, index=False, sep='\t')
motif_ls = list(motif_df[2])
site_num_ls2 = map(int, list(motif_df[4]))
pvalue_ls2 = [1.0e-100 if i < 1.0e-100 else i for i in map(float, list(motif_df[6]))]
tomtom_file = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/dynamic_region/egg_cell1/utr3/tomtom_egg_1cell_min5.tsv'

plot_motif_logo(motif_ls, site_num_ls2, pvalue_ls2, tomtom_file, rbp_name)