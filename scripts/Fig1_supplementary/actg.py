import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from pandas.io.parsers import read_csv
import seaborn as sns

a = read_csv('/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/ATCG/total.csv', index_col=0)
b = a[['A','T','C','G']]
b['total'] = b.apply(sum, axis=1)
c = b.apply(lambda x: x/x['total'] * 100, axis=1)
sample_ls = ['egg', '1cell', '4cell', '64cell', 'sphere', 'shield' ]
c = c.loc[sample_ls,['A','T','C','G']]
colors = [(138, 123, 176), (94, 149, 197), (160, 194, 121), (200, 107, 102)]
colors = ['#%02x%02x%02x'%i for i in colors]
sns.set_style('ticks')
c.columns = ['A','U','C','G']

c.index = ['0 h.p.f', '0.2 h.p.f', '1 h.p.f', '2 h.p.f', '4 h.p.f', '6 h.p.f']

draw_ls = ['0 h.p.f', '0.2 h.p.f', '1 h.p.f', '2 h.p.f', '4 h.p.f', '6 h.p.f']
fig, ax = plt.subplots(figsize=(3,3))
ax = c.loc[draw_ls,].plot(kind='bar',stacked=True, ylim=(0,100), color=colors, ax=ax )
ax.legend(bbox_to_anchor=(1.3, 1), borderaxespad=0., fontsize=11)
plt.ylabel('Percentage(%)', fontsize=14)
# plt.xlabel('Stages', fontsize=14)
plt.xticks(rotation=45, fontsize=12)
plt.yticks(fontsize=12)
ax.tick_params(top=False)
ax.tick_params(right=False)
# plt.subplots_adjust(right=0.9)
fig.savefig('actg.pdf', dpi=200, bbox_inches = "tight")