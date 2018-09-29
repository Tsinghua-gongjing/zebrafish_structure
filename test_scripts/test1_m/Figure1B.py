import pandas as pd
from pandas.io.parsers import read_csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def Figure1B():
    all_gini = '/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/global_str_change/data/RPKM_combine.merge.t200.mean_reactivity.null40.txt'
    df_all_gini = pd.read_csv(all_gini, header=0, index_col=0, sep='\t', na_values='NULL')
    sample_ls = ['egg', '1cell', '4cell', '64cell', 'sphere', 'shield']
    sample_labels = ['Egg', '1 Cell', '4 Cell', '64 Cell', 'Sphere', 'Shield']
    color_ls = sns.color_palette("Set1", n_colors=len(sample_ls), desat=0.8)
#     draw_ls = ['egg', '1cell', '4cell', '64cell', 'sphere', 'shield']
    with sns.axes_style("ticks"):
        fig, ax = plt.subplots(figsize=(6,4))
    feature = 'transcript'
    data = df_all_gini
    col_selected = ['%s(%s)'%(j,feature) for j in sample_ls]
    all_data = [list(data[i].dropna()) for i in col_selected]
    plot_data = pd.DataFrame()
    plot_df_ls = []
    for col, sample in zip(col_selected, sample_ls):
        tmp_data = pd.DataFrame()
        tmp_data['gini'] = data[col].dropna()
        tmp_data['sample'] = sample
        plot_df_ls.append(tmp_data)
    plot_data = pd.concat(plot_df_ls)
#     sns.pointplot(x='sample',y='gini', data=plot_data, estimator=np.median, ci=5, errwidth=1)
    x = [0, 0.5, 1, 2, 4, 6]
#     x = range(1, len(all_data)+1)
    medians = map(np.median, all_data)
    box = ax.boxplot(all_data, positions=x, labels=sample_labels, patch_artist=True, showbox=False, showfliers=False)
    for n,(patch,color) in enumerate(zip(box['boxes'], color_ls)):
        patch.set_facecolor(color)
        patch.set_edgecolor(color)
#     plt.setp(box['boxes'], visible=False)
    plt.setp(box['whiskers'], color='k', linestyle='-')
    plt.setp(box['medians'], visible=False)
    ax.set_ylim((0.1,0.5))
    ax.tick_params(top=False)
    ax.tick_params(right=False)
    sns.despine(top=True, right=True)
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        tick.set_fontsize(10)
    for tick in ax.get_yticklabels():
        tick.set_size(10)
    ax.plot(x, medians, linewidth=1, color='grey', linestyle='-')
    ax.scatter(x, medians, s=100, c=color_ls)
#     ax.set_xlabel('Development Stages', fontsize=14)
    ax.set_ylabel('mean icSHAPE Reactivity',fontsize=14)
    plt.tight_layout()
    plt.savefig("Figure2B.pdf", dpi=200, bbox_inches='tight')

if __name__ == '__main__':
    Figure1B()