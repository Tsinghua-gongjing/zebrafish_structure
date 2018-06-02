import sys,os
import pandas as pd
import numpy as np
from nested_dict import nested_dict
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")

def file_info(file_dir=None, result_dir=None):
    if file_dir is None:
        file_dir = '/Share/home/zhangqf7/gongjing/zebrafish/data/RT'
    if result_dir is None:
        result_dir = '/Share/home/zhangqf7/gongjing/zebrafish/result/RTCorrelationPairwise'

    files = os.listdir(file_dir)
    NAI_files = [i for i in files if i.startswith('NAI')]
    DMSO_files = [i for i in files if i.startswith('DMSO')]
    paths = [file_dir+'/'+i for i in files]
    file_info_dict = nested_dict()

    file_info_dict['file_dir'] = file_dir
    file_info_dict['files'] = files
    file_info_dict['paths'] = paths
    file_info_dict['result_dir'] = result_dir

    return file_info_dict.to_dict()

def cumulate_dist_plot(ls_ls, color_ls=None, ls_ls_label=None,bins=40,title=None,ax=None,savefn=None,xlabel=None,ylabel=None,add_vline=None,add_hline=None,log2transform=0,xlim=None,ylim=None):

    if log2transform:
        ls_ls = np.log2(ls_ls)
    values,base = np.histogram(ls_ls,bins=bins)
    cumulative = np.cumsum(values)
    cumulative_norm = [i/float(len(ls_ls)) for i in cumulative]
    if ls_ls_label:
        ax.plot(base[:-1],cumulative_norm,color=color_ls, label=ls_ls_label, linewidth=4)
    else:
        ax.plot(base[:-1],cumulative_norm,color=color_ls,linewidth=4)

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
            text.set_fontsize(20)
    plt.tight_layout()
    if savefn is not None:
        plt.savefig(savefn)
        plt.close()

def icshape_correlationRT_plot_combined_all_stages():
    file_info_dict = file_info(file_dir="/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/RT_corr_acc/data/",
                              result_dir="/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/RT_corr_acc/result/")
    sample_ls = ['egg', '1cell', '4cell', '64cell', 'sphere', 'shield', ]
    sample_labels = ['0 h.p.f', '0.2 h.p.f', '1 h.p.f', '2 h.p.f', '4 h.p.f', '6 h.p.f']
    sample_label_dict = dict(zip(sample_ls, sample_labels))
    color_ls = sns.color_palette("Set1", n_colors=len(sample_ls), desat=0.8)
    color_dict = dict(zip(sample_ls, color_ls))
    fig, ax = plt.subplots(2,6,sharex=True, sharey=True,figsize=(21,7))
    
    i = 0
    j = 0
    for group in ['DMSO', 'NAI']:
        print "group %s" % group
        for baseDensity_cutoff in [200]:
            print "plot cutoff %s" % baseDensity_cutoff
            outputDir = os.path.join(file_info_dict['result_dir'], "cutoff_T2t%s" % baseDensity_cutoff)
            if not os.path.exists(outputDir):
                os.makedirs(outputDir)
            j = 0
            for sample in sample_ls:
                sample_fn = group + "_" + sample + "_rep1" + "-" + group + "_" + sample + "_rep2" + "-" + "T2t%s" % baseDensity_cutoff
                if not os.path.exists(os.path.join(file_info_dict['file_dir'], sample_fn)):
                    print "file not exists: %s" % os.path.join(file_info_dict['file_dir'], sample_fn)
                    continue
                df = pd.read_csv(os.path.join(file_info_dict['file_dir'], sample_fn), header=None, sep='\t',  keep_default_na=False, na_values=['n/a'])
                df.dropna(axis=0, how='any',  inplace=True)
                print "sample read in %s " % sample_fn
                ls_ls = list(df[4])
                c_ls = color_dict.get(sample, 'black')

                ls_ls_label = "%s %s" % (sample_label_dict.get(sample), group)
                if not len(ls_ls) > 0:
                    continue

                
                cumulate_dist_plot(ax = ax[i,j], ls_ls=ls_ls, color_ls=c_ls, ls_ls_label=ls_ls_label, bins=60, title=None, savefn=None, xlabel=None, ylabel='', add_vline=None,add_hline=[0,1], log2transform=0, xlim=[-0.05,1.05], ylim=[-0.05,1.05])
                ax[i,j].tick_params(top=False)
                ax[i,j].tick_params(right=False)
                j = j + 1
#                 
        i = 1 + i
    plt.tight_layout()
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel("Pearson Correlation", fontsize=20)
    plt.ylabel("Percentage", fontsize=20)
    fig.savefig('rtcorr.pdf', dpi=200, bbox_inches = "tight")

def main():
	icshape_correlationRT_plot_combined_all_stages()

if __name__ == '__main__':
	main()