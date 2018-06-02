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
from scipy import stats

def read_rpkm_txt(txt, min_val=-1):
	gj.printFuncRun('read_rpkm_txt')
	gj.printFuncArgs()
	val_dict = nested_dict()
	gene_ls = []
	with open(txt, 'r') as TXT:
		for line in TXT:
			line = line.strip()
			if not line or line.startswith('#'): continue
			arr = line.split('\t')
			val_dict[arr[0]] = float(arr[4])
			gene_ls.append(arr[0])
	gj.printFuncRun('read_rpkm_txt')
	return val_dict,gene_ls

def corr_plot(fn_str=None, label_str=None, savefn=None, intersect='common'):
	gj.printFuncRun('corr_plot')
	gj.printFuncArgs()
	fn_ls = fn_str.split(':')
	label_ls = label_str.split(':')
	gene_ls_ls = []
	rpkm_dict_ls = []
	rpkm_ls_ls = []
	for fn in fn_ls:
		rpkm_dict,gene_ls = read_rpkm_txt(fn)
		rpkm_dict_ls.append(rpkm_dict)
		gene_ls_ls.append(gene_ls)
		rpkm_ls_ls.append([np.log2(i) for i in rpkm_dict.values()])

	gj.cumulate_dist_plot(ls_ls=rpkm_ls_ls,ls_ls_label=label_ls,bins=40,title=None,ax=None,savefn=savefn+'.cdf.png',xlabel='log2(RPKM)',ylabel=None,add_vline=None,add_hline=None,log2transform=0)

	if intersect == 'common':
		genes = gj.ls_ls_common(ls_ls=gene_ls_ls,return_ls=1)
	elif intersect == 'union':
		genes = gj.ls_ls_union(ls_ls=gene_ls_ls,return_ls=1)

	SAVEFN = open(savefn, 'w')
	print >>SAVEFN,'#gene'+'\t'+'\t'.join(label_ls)
	for gene in genes:
		gene_rpkm_ls = []
		for sample_rpkm_dict in rpkm_dict_ls:
			rpkm = np.log2(float(sample_rpkm_dict[gene])) if sample_rpkm_dict.has_key(gene) else np.log2(0.001)
			gene_rpkm_ls.append(rpkm)
		print >>SAVEFN,gene+'\t'+'\t'.join(map(str, gene_rpkm_ls))
	SAVEFN.close()

	df = pd.read_csv(savefn, sep='\t', header=0)
	gj.df_corr_matrix_plot(df[label_ls],savefn=savefn+'.png',size=4,rot=30,share_x_y=1,hue=None,diag='kde')

	r, p = stats.spearmanr(df['rep1'], df['rep2'])
	fig,ax = plt.subplots(figsize=(6,6))
	df.plot(kind='scatter', x='rep1', y='rep2', ax=ax)
	coor = 'spearman'
	ax.annotate("{} r = {:.2f}".format(coor,r),
                xy=(.3, .9),ha='center',va='center',fontsize=20, xycoords=ax.transAxes)
	plt.tight_layout()
	plt.savefig(savefn+'.2.png')
	plt.close()

	gj.printFuncRun('corr_plot')

if __name__ == "__main__":
	corr_plot(fn_str=sys.argv[1], label_str=sys.argv[2], savefn=sys.argv[3])
	"""
	python rpkm_corr.py sample1_label:sample2_label:sample3_label sample1_rpkm_file:sample2_rpkm_file:sample3_rpkm_file
	"""
