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

fig, ax = plt.subplots(figsize=(6,6))

size = 0.4

val1 = [9251937, 8605104, 7445960, 7840986] # transcriptome
val2 = [4310065, 4020865, 3461761, 3679722] # icSHAPE

cmap = plt.get_cmap("Set1")
outer_colors = cmap(np.arange(4)*4)
# inner_colors = cmap(np.array([1, 2, 5, 6, 9, 10]))

ax.pie(val1, radius=1, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w'), labels=['A', 'T', 'C', 'G'],autopct='%1.2f%%')

ax.pie(val2, radius=1-size,colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w'),autopct='%1.2f%%')

ax.set(aspect="equal")
plt.savefig('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/F1.all_structure_base_ratio.pie_same_radius.pdf')
plt.show()