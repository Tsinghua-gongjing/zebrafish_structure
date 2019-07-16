from pandas.io.parsers import read_csv
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

def Figure1D():
    all_gini = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/RPKM_combine.merge.t200.mean_reactivity.null40.txt'
    df_all_gini = pd.read_csv(all_gini, header=0, index_col=0, sep='\t', na_values='NULL')
    sample_ls = ['egg', '1cell', '4cell', '64cell', 'sphere', 'shield']
    sample_labels = ['Egg', '1 Cell', '4 Cell', '64 Cell', 'Sphere', 'Shield']
    color_ls = sns.color_palette("Set1", n_colors=len(sample_ls), desat=0.8)
    with sns.axes_style("ticks"):
        fig, ax = plt.subplots(figsize=(8,3))
    feature = 'transcript'
    data = df_all_gini
    col_selected = ['%s(%s)'%(j,feature) for j in sample_ls]
    data = data[col_selected].dropna(how='all')
    sample_pos = np.array([0, 0.5, 1, 2, 4, 6])
    data.columns = sample_labels
    data['tx_id'] = data.index
    # plot_data = pd.melt(data, id_vars=['tx_id'], value_vars=sample_pos, var_name='stage', value_name='gini')
    plot_data = [ data[i].dropna().values for i in sample_labels]
    x = [0, 0.5, 1, 2, 4, 6]
    # ax = sns.violinplot(x='stage', y='gini', data=plot_data, ax=ax, palette=color_ls, linewidth=1, cut=0,saturation=1, order=np.arange(0,6.2,0.2))
    violins = ax.violinplot(plot_data, positions=sample_pos, widths=0.22, showextrema=False)
    for n,(patch,color) in enumerate(zip(violins['bodies'], color_ls)):
            patch.set_facecolor(color)
            patch.set_edgecolor('k')
            patch.set_linewidth(1)
            patch.set_alpha(1)
    boxes = ax.boxplot(plot_data, positions=sample_pos,widths=0.03,patch_artist=True, showcaps=False, showfliers=False,labels=sample_labels)
    for n,(patch,color) in enumerate(zip(boxes['boxes'], color_ls)):
            patch.set_facecolor('k')
            patch.set_edgecolor('k')
            patch.set_alpha(1)
    plt.setp(boxes['medians'], visible=False)
    plt.setp(boxes['whiskers'], linewidth=0.5)
    medians = map(np.median, plot_data)
    ax.plot(x, medians, 'o', linewidth=1, color='grey', linestyle='-', markersize=3, markeredgecolor='grey', markerfacecolor='white', zorder=10)
    for xtk in ax.get_xticklabels():
        xtk.set_rotation(45)
        xtk.set_fontsize(10)
    # ax.scatter(x, medians, color='white', s=10)
    ax.set_ylabel('mean icshape reactivity')
    ax.set_ylim((0.2,0.4))
    ax.set_xlim((-0.2,6.2))
    sns.despine(top=True, right=True)
    print zip(sample_labels,medians)
    plt.savefig("/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/F1.all_sample_mean_reactivity_plot.pdf", dpi=200, bbox_inches='tight')

if __name__ == '__main__':
    Figure1D()
