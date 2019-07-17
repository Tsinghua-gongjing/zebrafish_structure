# -*- coding: utf-8 -*-

import gj
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
plt.rcParams["font.family"] = "helvetica"
import sys, os
from nested_dict import nested_dict
import numpy as np, pandas as pd
from scipy import stats

def plot(txt='/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/iCLIP/bam_corr/h4_full.readCounts.tab', sample1="'rep1'", sample2="'rep2'"):
	df = pd.read_csv(txt, header=0, sep='\t')
	df['log2(%s)'%(sample1)] = np.log2(df[sample1]+1)
	df['log2(%s)'%(sample2)] = np.log2(df[sample2]+1)
	print df.head()
	print df.columns

	fig,ax=plt.subplots(figsize=(8,8))
	# sns.jointplot(x=sample1, y=sample2, data=df, stat_func=stats.pearsonr, kind="reg")
	# sns.jointplot(x='log2(%s)'%(sample1), y='log2(%s)'%(sample2), data=df, stat_func=stats.pearsonr, kind="reg", size=10)
	df.plot(kind='scatter', x='log2(%s)'%(sample1), y='log2(%s)'%(sample2), s=2, ax=ax)
	r,p = stats.pearsonr(df['log2(%s)'%(sample1)], df['log2(%s)'%(sample2)])
	plt.title('r: %s, p: %s'%(r, p))
	lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
	]

	# now plot both limits against eachother
	ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
	ax.set_aspect('equal')
	ax.set_xlim(lims)
	ax.set_ylim(lims)

	savefn = (txt+'.%s.%s.joint.pdf'%(sample1, sample2)).replace("'", "")
	plt.savefig(savefn)
	plt.close()

plot()