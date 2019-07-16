from scipy import stats
import pandas as pd
from pandas.io.parsers import read_csv
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
plt.rcParams["font.family"] = "helvetica"


"""
/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail

head -1 /Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelationPairwiseNew/2018_03_12_gini/RPKM_combine.merge.addTwoRP33.t200.gini.null40.txt > 20190425_riboseq_DE.bin1000.upgenes.gini.txt
awk 'NR==FNR{a[$1];next} $1 in a{print $0}' 20190425_riboseq_DE.bin1000.upgenes.txt /Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelationPairwiseNew/2018_03_12_gini/RPKM_combine.merge.addTwoRP33.t200.gini.null40.txt >> 20190425_riboseq_DE.bin1000.upgenes.gini.txt

head -1 /Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelationPairwiseNew/2018_03_12_gini/RPKM_combine.merge.addTwoRP33.t200.gini.null40.txt > 20190425_riboseq_DE.bin1000.downgenes.gini.txt
awk 'NR==FNR{a[$1];next} $1 in a{print $0}' 20190425_riboseq_DE.bin1000.downgenes.txt /Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelationPairwiseNew/2018_03_12_gini/RPKM_combine.merge.addTwoRP33.t200.gini.null40.txt >> 20190425_riboseq_DE.bin1000.downgenes.gini.txt
"""

# def plot_dis(fn, sample1='1cell', sample2='1cell-RK-33'):
# 	df = pd.read_csv(fn, header=0, sep='\t')
# 	fig,ax=plt.subplots(1,4, sharex=False, sharey=False, figsize=(6*4, 6))
# 	df.plot(kind='scatter', x='%s(transcript)'%(sample1), y='%s(transcript)'%(sample2), ax=ax[0])
# 	df.plot(kind='scatter', x='%s(UTR5)'%(sample1), y='%s(UTR5)'%(sample2), ax=ax[1])
# 	df.plot(kind='scatter', x='%s(CDS)'%(sample1), y='%s(CDS)'%(sample2), ax=ax[2])
# 	df.plot(kind='scatter', x='%s(UTR3)'%(sample1), y='%s(UTR3)'%(sample2), ax=ax[3])

# 	for a in ax[0:]:
# 		lims = [
#     	np.min([a.get_xlim(), a.get_ylim()]),  # min of both axes
#     	np.max([a.get_xlim(), a.get_ylim()]),  # max of both axes
# 		]

# 		# now plot both limits against eachother
# 		a.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
# 		a.set_aspect('equal')
# 		a.set_xlim(lims)
# 		a.set_ylim(lims)

# 	plt.tight_layout()
# 	plt.savefig(fn.replace('.txt', '.element.pdf'))
# 	plt.close()

# plot_dis('/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/20190425_riboseq_DE.bin1000.downgenes.mean_reactivity.txt')

# def up_down():
# 	up = '/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/20190425_riboseq_DE.bin1000.upgenes.mean_reactivity.txt'
# 	down = '/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/20190425_riboseq_DE.bin1000.downgenes.mean_reactivity.txt'
# 	sample_ls = ['1cell(transcript)', '1cell-RK-33(transcript)', '1cell-RK-33-2(transcript)']
# 	sample_labels = ['1 Cell', '1 Cell (RK33-1)', '1 Cell (RK33-2)']
# 	df_up_gini = pd.read_csv(up, header=0, index_col=0, sep='\t', na_values='NULL')
# 	df_up_gini['type'] = 'up'
# 	df_down_gini = pd.read_csv(down, header=0, index_col=0, sep='\t', na_values='NULL')
# 	df_down_gini['type'] = 'down'
# 	# df_up_down_gini = df_up_gini.merge(df_down_gini, left_index=True, right_index=True, how='outer', suffixes=('up', 'down'))
# 	df_up_down_gini = pd.concat([df_up_gini, df_down_gini])[sample_ls+['type']]
# 	print df_up_gini.shape, df_down_gini.shape, df_up_down_gini.shape
# 	print df_up_down_gini.head()

# 	df_up_down_gini_melt = df_up_down_gini.melt(id_vars=['type'], value_vars=sample_ls, var_name='condition', value_name='gini')
# 	print df_up_down_gini_melt.head()

# 	fig,ax=plt.subplots(figsize=(8,12))
# 	g=sns.boxplot(x='condition', y='gini', hue='type', data=df_up_down_gini_melt)
# 	# g.set_xticklabels(rotation=45)
# 	ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=45)
# 	plt.tight_layout()
# 	plt.savefig('/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/up_down_mean_reactivity.png')
# 	plt.close()
# up_down()

# def total_down():
# 	total = '/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/total.RPKM_combine.merge.addTwoRP33.t200.gini.null40.txt'
# 	down = '/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/20190425_riboseq_DE.bin1000.downgenes.gini.txt'
# 	sample_ls = ['conrtol/rk33', 'conrtol/rk33-2']
# 	df_total_gini = pd.read_csv(total, header=0, index_col=0, sep='\t', na_values='NULL')
# 	df_total_gini['type'] = 'total'
# 	df_total_gini['conrtol/rk33'] = df_total_gini['1cell(UTR5)'] / df_total_gini['1cell-RK-33(UTR5)']
# 	df_total_gini['conrtol/rk33-2'] = df_total_gini['1cell(UTR5)'] / df_total_gini['1cell-RK-33-2(UTR5)']
# 	df_down_gini = pd.read_csv(down, header=0, index_col=0, sep='\t', na_values='NULL')
# 	df_down_gini['conrtol/rk33'] = df_down_gini['1cell(UTR5)'] / df_down_gini['1cell-RK-33(UTR5)']
# 	df_down_gini['conrtol/rk33-2'] = df_down_gini['1cell(UTR5)'] / df_down_gini['1cell-RK-33-2(UTR5)']
# 	df_down_gini['type'] = 'down'
# 	df_total_down_gini = pd.concat([df_total_gini, df_down_gini])[sample_ls+['type']]
# 	print df_total_gini.shape, df_down_gini.shape, df_total_down_gini.shape
# 	df_total_down_gini_melt = df_total_down_gini.melt(id_vars=['type'], value_vars=sample_ls, var_name='condition', value_name='gini(control)/gin(rk33)')
# 	print df_total_down_gini_melt.head()

# 	fig,ax=plt.subplots(figsize=(8,12))
# 	g=sns.boxplot(x='condition', y='gini(control)/gin(rk33)', hue='type', data=df_total_down_gini_melt)
# 	# g.set_xticklabels(rotation=45)
# 	ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=45)
# 	plt.tight_layout()
# 	plt.savefig('/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/total_down_mean_reactivity.png')
# 	plt.close()

# total_down()

def structure_diff_vs_TE_diff(structure=None, TE=None, label1='control', label2='RK33', feature='UTR5'):
	if structure is None:
		structure = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/RPKM_combine.merge.addTwoRP33.t200.gini.null40.txt'
	if TE is None:
		TE = '/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/20190425_TE_corr.txt'

	df_TE = pd.read_csv(TE, header=0, sep='\t')
	df_TE['diff(TE, minus)'] = df_TE['mean(TE(%s))'%(label1)] - df_TE['mean(TE(%s))'%(label2)] 
	df_TE['log2[TE(%s)-TE(%s)]'%(label1, label2)] = np.log2(df_TE['diff(TE, minus)'])
	df_TE['diff(TE, divide)'] = df_TE['mean(TE(%s))'%(label1)] / df_TE['mean(TE(%s))'%(label2)] 
	df_TE['log2[TE(%s)/TE(%s)]'%(label1, label2)] = np.log2(df_TE['diff(TE, divide)'])

	df_structure = pd.read_csv(structure, header=0, sep='\t')
	df_structure['control(%s)'%(feature)] = df_structure['1cell(%s)'%(feature)]
	df_structure['RK33(%s)'%(feature)] = df_structure['1cell-RK-33(%s)'%(feature)]
	df_structure['transcript'] = df_structure.index
	df_structure['Gini(%s)-Gini(%s)'%(label1, label2)] = df_structure['%s(%s)'%(label1, feature)] - df_structure['%s(%s)'%(label2, feature)]
	df_structure['Gini(%s)/Gini(%s)'%(label1, label2)] = df_structure['%s(%s)'%(label1, feature)] / df_structure['%s(%s)'%(label2, feature)]
	# df_structure['conrtol/rk33-2'] = df_structure['1cell(transcript)'] - df_structure['1cell-RK-33-2(transcript)']

	df_structure_TE = df_structure.merge(df_TE, on='transcript', how='inner')
	# df_structure_TE = df_structure_TE[df_structure_TE['log2(diff(TE))']<10]
	print df_structure_TE.shape
	print df_structure_TE.columns

	fig,ax=plt.subplots(figsize=(10,10))
	g=sns.jointplot(x='Gini(%s)/Gini(%s)'%(label1, label2), y='log2[TE(%s)-TE(%s)]'%(label1, label2), data=df_structure_TE, kind="reg", stat_func=stats.pearsonr, size=10)
	n=df_structure_TE[['Gini(%s)/Gini(%s)'%(label1, label2), 'log2[TE(%s)-TE(%s)]'%(label1, label2)]].dropna(how='any').shape[0]
	a = g.fig.get_axes()
	a[2].set_title('n=%s'%(n))
	plt.tight_layout()
	plt.savefig('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/riboseq/structure_diff_vs_TE_diff.GiniDivideTEminus.pdf')
	plt.close()

	# fig,ax=plt.subplots(figsize=(10,10))
	# g=sns.jointplot(x='Gini(%s)/Gini(%s)'%(label1, label2), y='log2[TE(%s)/TE(%s)]'%(label1, label2), data=df_structure_TE, kind="reg", stat_func=stats.pearsonr, size=10)
	# n=df_structure_TE[['Gini(%s)/Gini(%s)'%(label1, label2), 'log2[TE(%s)/TE(%s)]'%(label1, label2)]].dropna(how='any').shape[0]
	# a = g.fig.get_axes()
	# a[2].set_title('n=%s'%(n))
	# plt.tight_layout()
	# plt.savefig('/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/structure_diff_vs_TE_diff.GiniDivideTEDivide.pdf')
	# plt.close()

	# fig,ax=plt.subplots(figsize=(10,10))
	# g=sns.jointplot(x='Gini(%s)-Gini(%s)'%(label1, label2), y='log2[TE(%s)-TE(%s)]'%(label1, label2), data=df_structure_TE, kind="reg", stat_func=stats.pearsonr, size=10)
	# n=df_structure_TE[['Gini(%s)-Gini(%s)'%(label1, label2), 'log2[TE(%s)-TE(%s)]'%(label1, label2)]].dropna(how='any').shape[0]
	# a = g.fig.get_axes()
	# a[2].set_title('n=%s'%(n))
	# plt.tight_layout()
	# plt.savefig('/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/structure_diff_vs_TE_diff.GiniMinusTEminus.pdf')
	# plt.close()

	# fig,ax=plt.subplots(figsize=(10,10))
	# g=sns.jointplot(x='Gini(%s)-Gini(%s)'%(label1, label2), y='log2[TE(%s)/TE(%s)]'%(label1, label2), data=df_structure_TE, kind="reg", stat_func=stats.pearsonr, size=10)
	# n=df_structure_TE[['Gini(%s)-Gini(%s)'%(label1, label2), 'log2[TE(%s)/TE(%s)]'%(label1, label2)]].dropna(how='any').shape[0]
	# a = g.fig.get_axes()
	# a[2].set_title('n=%s'%(n))
	# plt.tight_layout()
	# plt.savefig('/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/structure_diff_vs_TE_diff.GiniMinusTEDivide.pdf')
	# plt.close()
	

	# for i in ['transcript', 'UTR5', 'CDS', 'UTR3']:
	# 	fig,ax=plt.subplots()
	# 	g=sns.jointplot(x='mean(TE(control))', y='1cell(%s)'%(i), data=df_structure_TE, kind="reg", stat_func=stats.pearsonr, size=10)
	# 	plt.tight_layout()
	# 	plt.savefig('/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/structure_vs_TE_1cell.%s.png'%(i))
	# 	plt.close()
	# for i in ['transcript', 'UTR5', 'CDS', 'UTR3']:
	# 	fig,ax=plt.subplots()
	# 	g=sns.jointplot(x='mean(TE(RK33))', y='1cell-RK-33(%s)'%(i), data=df_structure_TE, kind="reg", stat_func=stats.pearsonr, size=10)
	# 	plt.tight_layout()
	# 	plt.savefig('/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/xtail/structure_vs_TE_1cell_RK33.%s.png'%(i))
	# 	plt.close()

structure_diff_vs_TE_diff(TE='/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/riboseq/20190415_TE_corr.txt', feature='transcript')
# structure_diff_vs_TE_diff(structure=None, TE='/Share/home/zhangqf7/gongjing/zebrafish/result/riboseq_rnaseq/20190425_TE_egg_1cell_corr.txt', label1='egg', label2='1cell', feature='UTR5')
# montage structure_diff_vs_TE_diff.G*png  -mode concatenate -tile 2x structure_diff_vs_TE_diff.pdf


