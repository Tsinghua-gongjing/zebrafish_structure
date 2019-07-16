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
import gj

def TE_rep_corr(TE1, TE2, savefn, label1='control', label2='RK33'):
	df_TE1 = pd.read_csv(TE1, header=0, sep='\t')
	df_TE2 = pd.read_csv(TE2, header=0, sep='\t')

	fig,ax=plt.subplots()
	ls_ls = [list(df_TE2['log2(TE(%s))'%(label1)]), list(df_TE2['log2(TE(%s))'%(label2)])]
	# df_TE1[['log2(TE(control))', 'log2(TE(RK33))']].plot(kind='bar', ax=ax)
	gj.cumulate_dist_plot(ls_ls=ls_ls,ls_ls_label=[label1, label2],bins=40000,title=None,ax=None,savefn=TE1+'.cumulate.png',xlabel=None,ylabel=None,add_vline=None,add_hline=None,log2transform=0,xlim=None,ylim=None)

	df_merge = df_TE1.merge(df_TE2, on='transcript', how='inner')

	df_merge['mean(TE(%s))'%(label1)] = (df_merge['TE(%s)_x'%(label1)] + df_merge['TE(%s)_y'%(label1)]) / 2.0
	df_merge['mean(TE(%s))'%(label2)] = (df_merge['TE(%s)_x'%(label2)] + df_merge['TE(%s)_y'%(label2)]) / 2.0
	df_merge['log2(mean(TE(%s)))'%(label1)] = np.log2(df_merge['mean(TE(%s))'%(label1)])
	df_merge['log2(mean(TE(%s)))'%(label2)] = np.log2(df_merge['mean(TE(%s))'%(label2)])
	ls_ls = [list(df_merge['log2(mean(TE(%s)))'%(label1)]), list(df_merge['log2(mean(TE(%s)))'%(label2)])]
	p = gj.ks_2samp(ls_ls[0], ls_ls[1])
	print "pvalue: %s"%(p)
	gj.cumulate_dist_plot(ls_ls=ls_ls,ls_ls_label=[label1, label2],bins=40000,title=None,ax=None,savefn=TE1+'.cumulate.mean.pdf',xlabel=None,ylabel=None,add_vline=None,add_hline=None,log2transform=0,xlim=[-5, 5],ylim=None)

	df_merge.to_csv(savefn, header=True, index=False, sep='\t')

	fig,ax = plt.subplots(figsize=(8,8))
	df_merge.plot(kind='scatter', x='log2(TE(%s))_x'%(label1), y='log2(TE(%s))_y'%(label1), ax=ax)
	r,p = stats.pearsonr(df_merge['log2(TE(%s))_x'%(label1)], df_merge['log2(TE(%s))_y'%(label1)])
	plt.title("r: %s, p:%s"%(r, p))

	lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
	]
	
	# now plot both limits against eachother
	ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
	ax.set_aspect('equal')
	ax.set_xlim(lims)
	ax.set_ylim(lims)

	plt.tight_layout()
	plt.savefig(savefn.replace('.txt', '.%s.pdf'%(label1)))
	plt.close()

	fig,ax = plt.subplots(figsize=(8,8))
	df_merge.plot(kind='scatter', x='log2(TE(%s))_x'%(label2), y='log2(TE(%s))_y'%(label2), ax=ax)
	r,p = stats.pearsonr(df_merge['log2(TE(%s))_x'%(label2)], df_merge['log2(TE(%s))_y'%(label2)])
	plt.title("r: %s, p:%s"%(r, p))

	lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
	]
	
	# now plot both limits against eachother
	ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
	ax.set_aspect('equal')
	ax.set_xlim(lims)
	ax.set_ylim(lims)

	plt.tight_layout()
	plt.savefig(savefn.replace('.txt', '.%s.pdf'%(label2)))
	plt.close()


def main():
	TE1 = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/riboseq/20190415_riboseq_vs_rnaseq_control_vs_rk33_rep1.TE.txt'
	TE2 = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/riboseq/20190415_riboseq_vs_rnaseq_control_vs_rk33_rep2.TE.txt'
	savefn = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/riboseq/20190415_TE_corr.txt'
	TE_rep_corr(TE1, TE2, savefn)#, label1='RK33-1', label2='RK33-2')

if __name__ == '__main__':
	main()