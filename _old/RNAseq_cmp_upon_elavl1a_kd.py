import pandas as pd
import numpy as np
import os
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats


def cumulate_dist_plot(ls_ls, color_ls=None, ls_ls_label=None, bins=1000, ls='-', title=None, ax=None, savefn=None, xlabel=None, ylabel=None, add_vline=None, add_hline=None, log2transform=0, xlim=None, ylim=None):
    if ax is None:
        sns.set_style("ticks")
        fig, ax = plt.subplots(figsize=(6, 6))

    if log2transform:
        ls_ls = np.log2(ls_ls)
    values, base = np.histogram(ls_ls, bins=bins)
#     print values, base
    cumulative = np.cumsum(values)
    cumulative_norm = [i / float(len(ls_ls)) for i in cumulative]
    if ls_ls_label:
        ax.plot(base[:-1], cumulative_norm, color=color_ls,
                label=ls_ls_label, linewidth=4, ls=ls)
    else:
        ax.plot(base[:-1], cumulative_norm, color=color_ls, linewidth=4)
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
        ax.set_title(title, y=1.05)
    if ls_ls_label:
        pass

    if add_vline is not None:
        for vline in add_vline:
            ax.axvline(vline, ls="--", color='lightgrey')
    if add_hline is not None:
        for hline in add_hline:
            ax.axhline(hline, ls="--", color='lightgrey')
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


ctrl_6h_vs_4h_file = './RNAseq_ctrl_h46/output_score.txt'
ctrl_6h_vs_4h_df = pd.read_csv(ctrl_6h_vs_4h_file, header=0, sep="\t").dropna()
ctrl_6h_vs_4h_df.index = ctrl_6h_vs_4h_df['GeneNames']
mo_6h_vs_4h_file = './RNAseq-hur-mo-h46/output_score.txt'
mo_6h_vs_4h_df = pd.read_csv(mo_6h_vs_4h_file, header=0, sep="\t").dropna()
mo_6h_vs_4h_df.index = mo_6h_vs_4h_df['GeneNames']

decay = pd.read_csv('./decay.txt', sep="\t")
decay['tx_id'] = decay.index
decay['type'] = 'decay'
stable = pd.read_csv('./stable.txt', sep="\t")
stable['tx_id'] = stable.index
stable['type'] = 'stable'
zygotic = pd.read_csv('./zygotic.txt', sep="\t")
zygotic['tx_id'] = zygotic.index
zygotic['type'] = 'zygotic'
group = pd.concat([decay[['tx_id', 'type']],
                   stable[['tx_id', 'type']], zygotic[['tx_id', 'type']]])
ctrl_data = pd.concat([ctrl_6h_vs_4h_df, group], axis=1, join='outer').dropna()
mo_data = pd.concat([mo_6h_vs_4h_df, group], axis=1, join='outer').dropna()

ctrl_stable_group = ctrl_data[ctrl_data['type'] == 'stable']
ctrl_decay_group = ctrl_data[ctrl_data['type'] == 'decay']
ctrl_zygotic_group = ctrl_data[ctrl_data['type'] == 'zygotic']
mo_stable_group = mo_data[mo_data['type'] == 'stable']
mo_decay_group = mo_data[mo_data['type'] == 'decay']
mo_zygotic_group = mo_data[mo_data['type'] == 'zygotic']

sns.set_style('ticks')
fig, ax = plt.subplots()
stable_color = '#00008b'
decay_color = '#a52a2a'
zygotic_color = '#ffb202'
cumulate_dist_plot(ctrl_decay_group['log2(Fold_change) normalized'], color_ls=decay_color,
                   ls_ls_label='ctrl_decay (%d)' % len(ctrl_decay_group), ax=ax, ls='--')
cumulate_dist_plot(ctrl_stable_group['log2(Fold_change) normalized'], color_ls=stable_color,
                   ls_ls_label='ctrl_stable (%d)' % len(ctrl_stable_group), ax=ax, ls='--')
cumulate_dist_plot(ctrl_zygotic_group['log2(Fold_change) normalized'], color_ls=zygotic_color,
                   ls_ls_label='ctrl_zygotic (%d)' % len(ctrl_zygotic_group), ax=ax, ls='--')
cumulate_dist_plot(mo_decay_group['log2(Fold_change) normalized'],
                   color_ls=decay_color, ls_ls_label='mo_decay (%d)' % len(mo_decay_group), ax=ax)
cumulate_dist_plot(mo_stable_group['log2(Fold_change) normalized'],
                   color_ls=stable_color, ls_ls_label='mo_stable (%d)' % len(mo_stable_group), ax=ax)

cumulate_dist_plot(mo_zygotic_group['log2(Fold_change) normalized'], color_ls=zygotic_color,
                   ls_ls_label='mo_zygotic (%d)' % len(mo_zygotic_group), ax=ax, add_hline=[0, 1])
plt.xlim((-7, 7))
plt.xlabel('log(fc of rpkm, 6_vs_4)')
# print "ctrl_decay: %d" % len(ctrl_decay_group)
# print "ctrl_stable: %d" % len(ctrl_stable_group)
# print "ctrl_zygotic: %d" % len()
# print "mo_decay: %d" % len(mo_decay_group)
# print "mo_stable: %d" % len(mo_stable_group)
print "decay", stats.ks_2samp(ctrl_decay_group['log2(Fold_change) normalized'], mo_decay_group['log2(Fold_change) normalized'])
print "stable", stats.ks_2samp(ctrl_stable_group['log2(Fold_change) normalized'], mo_stable_group['log2(Fold_change) normalized'])
print "zygotic", stats.ks_2samp(ctrl_zygotic_group['log2(Fold_change) normalized'], mo_zygotic_group['log2(Fold_change) normalized'])
_, d_p = stats.ks_2samp(
    ctrl_decay_group['log2(Fold_change) normalized'], mo_decay_group['log2(Fold_change) normalized'])
_, s_p = stats.ks_2samp(
    ctrl_stable_group['log2(Fold_change) normalized'], mo_stable_group['log2(Fold_change) normalized'])
_, z_p = stats.ks_2samp(
    ctrl_zygotic_group['log2(Fold_change) normalized'], mo_zygotic_group['log2(Fold_change) normalized'])
ax.text(0.8, 0.2, 'p < %.2e' % (max(
    [d_p, s_p, z_p])), fontsize=12, ha='center', va='center', transform=ax.transAxes)
plt.savefig('./hur-mo_vs_wt_3g_decay_cum_plot.pdf',
            dpi=220, bbox_inches='tight')
