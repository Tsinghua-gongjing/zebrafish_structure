import sys,os
import pandas as pd
from pandas.io.parsers import read_csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")

df = read_csv('/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/transcript_num_with_str/20180526_stats.txt', sep="\t", header=None, index_col=0)
df.index = ['Egg', '1 Cell', '4 Cell', '64 Cell', 'Sphere', 'Shield']
sample_ls = ['Egg', '1 Cell', '4 Cell', '64 Cell', 'Sphere', 'Shield']
df_plot = df.loc[sample_ls,1]
color_ls = sns.color_palette("Set1", n_colors=len(df_plot), desat=0.8)
fig, ax = plt.subplots(figsize=(8,4))
sns.barplot(x=list(df_plot.index), y=df_plot.values, palette=color_ls, ax=ax)
ax.tick_params(top=False)
ax.tick_params(right=False)
sns.despine(top=True, right=True)
ax.set_ylabel("Transcript Number", fontsize=14)
plt.xticks(rotation=45, fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('transcript_number.pdf', dpi=200, bbox_inches = "tight")