# -*- coding: utf-8 -*-

import gj
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
plt.rcParams["font.family"] = "helvetica"
import sys
import pandas as pd 
import venn_all

fn_ls = sys.argv[1].split(':')
labels_ls = sys.argv[2].split(':')
savefn = sys.argv[3]

id_ls = []
for i in fn_ls:
	df = pd.read_csv(i, header=None, sep='\t')
	id_ls.append(set(df[0]))

# matplot venn: cannot set same equal size
# gj.venn3plot(mode='string', subsets_ls=id_ls, labels_ls=labels_ls, title_str='', save_fn=savefn)

labels = venn_all.get_labels(data=id_ls)
if len(labels_ls) == 2:
	venn_all.venn2(labels=labels, names=labels_ls)
if len(labels_ls) == 3:
	venn_all.venn3(labels=labels, names=labels_ls)
if len(labels_ls) == 4:
	venn_all.venn4(labels=labels, names=labels_ls, figsize=(6,6))
plt.tight_layout()
plt.savefig(savefn)
plt.close()
