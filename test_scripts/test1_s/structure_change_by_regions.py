import sys,os
import pandas as pd
from pandas.io.parsers import read_csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


all_gini = '/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/global_str_change/data/RPKM_combine.merge.t200.gini.null40.txt'
df_all_gini = pd.read_csv(all_gini, header=0, index_col=0, sep='\t', na_values='NULL')
all_m_ic = '/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/global_str_change/data/RPKM_combine.merge.t200.mean_reactivity.null40.txt'
df_all_m_ic = pd.read_csv(all_m_ic, header=0, index_col=0, sep='\t', na_values='NULL')
sample_ls = ['egg', '1cell', '4cell', '64cell', 'sphere', 'shield']
sample_labels = ['Egg', '1 Cell', '4 Cell', '64 Cell', 'Sphere', 'Shield']
color_ls = sns.color_palette("Set1", n_colors=len(sample_ls), desat=0.8)
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(2,3, sharey='row', sharex=True, figsize=(7,5))
k = 0
for data in [df_all_gini, df_all_m_ic]:
    l = 0
    for feature in ['UTR5', 'CDS', 'UTR3']:
        col_selected = ['%s(%s)'%(j,feature) for j in sample_ls]
        all_data = [list(data[i].dropna()) for i in col_selected]
        medians = map(np.median, all_data)
        x = range(len(all_data))
        box = ax[k,l].boxplot(all_data, positions=x, labels=sample_labels, patch_artist=True, showfliers=False )
        for n,(patch,color) in enumerate(zip(box['boxes'], color_ls)):
            patch.set_facecolor(color)
            patch.set_edgecolor(color)
    #     for.
        plt.setp(box['whiskers'], color='k', linestyle='-')
        plt.setp(box['medians'], color='k')
        if k == 0:
            ax[k,l].set_ylim((0.2,0.8)) # gini
        else:
            ax[k,l].set_ylim((0.1,0.7)) # mean ic
        ax[k,l].tick_params(top=False)
        ax[k,l].tick_params(right=False)
        for tick in ax[k,l].get_xticklabels():
            tick.set_rotation(45)
            tick.set_fontsize(10)
        for tick in ax[k,l].get_yticklabels():
            tick.set_size(10)
        ax[k,l].plot(x, medians, linewidth=1, color='grey', linestyle='-')
        l = l + 1
    k = k + 1
ax[0,0].set_ylabel('gini', fontsize=14)
ax[1,0].set_ylabel('Mean Reactivity', fontsize=14)
plt.tight_layout()
plt.savefig('str_change_5utrcds3utr.pdf', dpi=200)