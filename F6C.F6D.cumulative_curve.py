import pandas as pd
import numpy as np
import os
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats

def cumulate_dist_plot(ls_ls, color_ls=None, ls_ls_label=None,bins=1000,ls='-',title=None,ax=None,savefn=None,xlabel=None,ylabel=None,add_vline=None,add_hline=None,log2transform=0,xlim=None,ylim=None):
    if ax is None:
        sns.set_style("ticks")
        fig,ax = plt.subplots(figsize=(6,6))  
                                  
    if log2transform:
        ls_ls = np.log2(ls_ls)
    values,base = np.histogram(ls_ls,bins=bins)
#     print values, base
    cumulative = np.cumsum(values)
    cumulative_norm = [i/float(len(ls_ls)) for i in cumulative]
    from scipy.interpolate import spline
    T = np.array(base[:-1])
    xnew = np.linspace(T.min(),T.max(),3000) #300 represents number of points to make between T.min and T.max

    power_smooth = spline(T,cumulative_norm,xnew)

    ax.plot(xnew,power_smooth, color=color_ls,label='smooth_'+ls_ls_label+"(%d)" % len(ls_ls), ls=ls )
#     plt.show()
#     if ls_ls_label:
#         ax.plot(base[:-1],cumulative_norm,color=color_ls, label=ls_ls_label, linewidth=2, ls=ls)
#     else:
#         ax.plot(base[:-1],cumulative_norm,color=color_ls,linewidth=2)
    # data = np.array(data)
    # sorted_data = np.sort(data)
    # # sorted_data = sorted_data/ max(sorted_data)
    # ax.step(np.concatenate([sorted_data, sorted_data[[-1]]]), np.arange(sorted_data.size+1) / float(sorted_data.size+1), color=color_ls[n], label=ls_ls_label[n])

#     print "plot line num: %s"%(n)
    if ls_ls_label:    
        ax.legend(loc='best')
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel("Accumulate percent over total")
    if title is not None:
        ax.set_title(title, y = 1.05)
    if ls_ls_label:
        pass

    if add_vline is not None:
        for vline in add_vline:
            ax.axvline(vline,ls="--", color='lightgrey')
    if add_hline is not None:
        for hline in add_hline:
            ax.axhline(hline,ls="--", color='lightgrey')
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    if ls_ls_label:
        leg = ax.legend(loc="upper left")
        for text in leg.get_texts():
            text.set_fontsize(10)
    plt.tight_layout()
    if savefn is not None:
        plt.savefig(savefn)
        plt.close()

wd = '/Users/soul/BaiduyunDisk/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/FIGURE3/RNAseq_old/'
os.chdir(wd)

ctrl_6h_vs_4h_file = './RNAseq_ctrl_h46/output_score.txt'
ctrl_6h_vs_4h_df = pd.read_csv(ctrl_6h_vs_4h_file, header=0, sep="\t").dropna()
ctrl_6h_vs_4h_df.index = ctrl_6h_vs_4h_df['GeneNames']
mo_6h_vs_4h_file = './RNAseq-hur-mo-h46/output_score.txt'
mo_6h_vs_4h_df = pd.read_csv(mo_6h_vs_4h_file, header=0, sep="\t").dropna()
mo_6h_vs_4h_df.index = mo_6h_vs_4h_df['GeneNames']

decay = pd.read_csv('./decay.txt', sep="\t")
decay['tx_id'] = decay.index
decay['type'] = 'decay'
stable = pd.read_csv('./stable.txt',sep="\t")
stable['tx_id'] = stable.index
stable['type'] = 'stable'
zygotic = pd.read_csv('./zygotic.txt', sep="\t")
zygotic['tx_id'] = zygotic.index
zygotic['type'] = 'zygotic'
group = pd.concat([decay[['tx_id', 'type']], stable[['tx_id', 'type']], zygotic[['tx_id', 'type']]])
ctrl_data = pd.concat([ctrl_6h_vs_4h_df, group], axis=1, join='outer').dropna()
mo_data = pd.concat([mo_6h_vs_4h_df, group],axis=1, join='outer').dropna()

ctrl_stable_group = ctrl_data[ctrl_data['type'] == 'stable']
ctrl_decay_group = ctrl_data[ctrl_data['type'] == 'decay']
ctrl_zygotic_group = ctrl_data[ctrl_data['type'] == 'zygotic']
mo_stable_group = mo_data[mo_data['type'] == 'stable']
mo_decay_group = mo_data[mo_data['type'] == 'decay']
mo_zygotic_group = mo_data[mo_data['type'] == 'zygotic']

fig, ax = plt.subplots()
cumulate_dist_plot(ctrl_6h_vs_4h_df['log2(Fold_change) normalized'], bins=1000,color_ls='k', ls_ls_label='ctrl_all', ax=ax, ls='--')
cumulate_dist_plot(mo_6h_vs_4h_df['log2(Fold_change) normalized'], bins=1000,color_ls='grey', ls_ls_label='mo_all', ax=ax)
plt.xlim((-8,8))
plt.xlabel('log(fc of rpkm, 6_vs_4)')
print "all", stats.ks_2samp(ctrl_6h_vs_4h_df['log2(Fold_change) normalized'], mo_6h_vs_4h_df['log2(Fold_change) normalized'])
# print "all", stats.mannwhitneyu(ctrl_6h_vs_4h_df['log2(Fold_change) normalized'], mo_6h_vs_4h_df['log2(Fold_change) normalized'], alternative='greater')
print len(ctrl_6h_vs_4h_df['log2(Fold_change) normalized'])
print len(mo_6h_vs_4h_df['log2(Fold_change) normalized'])


##
shield_ctrl_vs_mo_file = './RNAseq-hur-mo-h6/output_score.txt'
shield_ctrl_vs_mo_df = pd.read_csv(shield_ctrl_vs_mo_file, header=0, sep="\t").dropna()
shield_ctrl_vs_mo_df.index = shield_ctrl_vs_mo_df['GeneNames']
shield_ctrl_vs_mo_data = pd.concat([shield_ctrl_vs_mo_df, group],axis=1,join='outer').dropna()
shield_ctrl_vs_mo_decay = shield_ctrl_vs_mo_data[shield_ctrl_vs_mo_data['type']=='decay'].dropna()
shield_ctrl_vs_mo_stable = shield_ctrl_vs_mo_data[shield_ctrl_vs_mo_data['type']=='stable'].dropna()
shield_ctrl_vs_mo_zygotic = shield_ctrl_vs_mo_data[shield_ctrl_vs_mo_data['type']=='zygotic'].dropna()
fig, ax = plt.subplots()
# cumulate_dist_plot(shield_ctrl_vs_mo_data['log2(Fold_change) normalized'],ax=ax,
#                   color_ls='r', ls_ls_label ='decay')
cumulate_dist_plot(shield_ctrl_vs_mo_stable['log2(Fold_change) normalized'],
                   color_ls='b', ls_ls_label='stable',bins=1000,ax=ax)
cumulate_dist_plot(shield_ctrl_vs_mo_decay['log2(Fold_change) normalized'],
                  color_ls='r', ls_ls_label='decay',bins=1000,ax=ax)
cumulate_dist_plot(shield_ctrl_vs_mo_zygotic['log2(Fold_change) normalized'],
                  color_ls='orange', ls_ls_label='zygotic', bins=1000, ax=ax)
plt.xlim((-2,2))
print shield_ctrl_vs_mo_stable['log2(Fold_change) normalized'].shape
print shield_ctrl_vs_mo_decay['log2(Fold_change) normalized'].shape
print len(np.unique(shield_ctrl_vs_mo_stable['log2(Fold_change) normalized']))
print len(np.unique(shield_ctrl_vs_mo_decay['log2(Fold_change) normalized']))
print stats.ks_2samp(shield_ctrl_vs_mo_stable['log2(Fold_change) normalized'], shield_ctrl_vs_mo_decay['log2(Fold_change) normalized'],)
# print stats.mannwhitneyu(shield_ctrl_vs_mo_stable['log2(Fold_change) normalized'], shield_ctrl_vs_mo_decay['log2(Fold_change) normalized'], 
                         # alternative='less')
plt.xlabel('6 h.p.f. mo/ctrl')