import sys,os
import pandas as pd
from pandas.io.parsers import read_csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def gini_plot():
    all_gini = '/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/global_str_change/data/RPKM_combine.merge.t200.gini.null40.txt'
    df_all_gini = pd.read_csv(all_gini, header=0, index_col=0, sep='\t', na_values='NULL')
    sample_ls = ['egg', '1cell', '4cell', '64cell', 'sphere', 'shield']
    sample_labels = ['Egg', '1 Cell', '4 Cell', '64 Cell', 'Sphere', 'Shield']
    color_ls = sns.color_palette("Set1", n_colors=len(sample_ls), desat=0.8)
    with sns.axes_style("ticks"):
        fig, ax = plt.subplots(figsize=(6,3))
    feature = 'transcript'
    data = df_all_gini
    col_selected = ['%s(%s)'%(j,feature) for j in sample_ls]
    all_data = [list(data[i].dropna()) for i in col_selected]
    x = [0, 0.5, 1, 2, 4, 6]
    medians = map(np.median, all_data)
    box = ax.boxplot(all_data, positions=x, labels=sample_labels, patch_artist=True,showfliers=False)
    for n,(patch,color) in enumerate(zip(box['boxes'], color_ls)):
        patch.set_facecolor(color)
        patch.set_edgecolor(color)
    plt.setp(box['boxes'], visible=False)
    plt.setp(box['whiskers'], color='k', linestyle='-')
    plt.setp(box['medians'], visible=False)
    ax.set_ylim((0.3,0.7))
    ax.tick_params(top=False)
    ax.tick_params(right=False)
    sns.despine(top=True, right=True)
    for tick in ax.get_xticklabels():
        tick.set_fontsize(10)
    for tick in ax.get_yticklabels():
        tick.set_size(10)
    ax.plot(x, medians, linewidth=1, color='grey', linestyle='-')
    ax.scatter(x, medians, s=100, c=color_ls)
    ax.set_xlabel('Development Stages', fontsize=14)
    ax.set_ylabel('gini index',fontsize=14)
    plt.tight_layout()
    plt.savefig('gini_all.pdf', dpi=200, )

if __name__ == '__main__':
    gini_plot()