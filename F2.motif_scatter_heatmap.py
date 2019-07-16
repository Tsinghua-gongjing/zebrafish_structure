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
from adjustText import adjust_text

f = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/dynamic_region/dynamic_region.way2345.UTR3.fimo.enrich.txt'
df = pd.read_csv(f, header=0, sep='\t')
df.head()

df['pvalue_adj'] = [10e-300 if float(i)==0 else i for i in df['pvalue_adj']]
df = df[df['pvalue_adj']<=0.05]
df['log2(fimo1_motif)'] = np.log2(df['fimo1_motif'])
df['-log10(pvalue_adj)'] = -np.log10(df['pvalue_adj'])

# print df

fig,ax = plt.subplots(figsize=(25,20))
df.plot(kind='scatter', x='log2(fimo1_motif)', y='-log10(pvalue_adj)', ax=ax)
ax.set_xlim(0,)

texts = []
for x,y,t in zip(df['log2(fimo1_motif)'], df['-log10(pvalue_adj)'], df['RBP_name']):
#     print x,y,t
    texts.append(plt.text(x, y, t, fontsize=20))
adjust_text(texts, only_move={'text': 'xy'})

plt.tight_layout()
plt.savefig('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/F2.dynamic_region.way2345.UTR3.fimo.enrich.pdf')


df_h20 = df.sort_values(by='pvalue_adj', ascending=True).head(100)
df_h20 = df_h20[df_h20['odd']>1]
df_h20 = df_h20[df_h20['RBP_name']!='HuR, ELAVL2/3/4']
df_h20['-log10(pvalue_adj)'] = -np.log10(df_h20['pvalue_adj'])
df_h20 = df_h20[df_h20['-log10(pvalue_adj)']>5]
df_h20.index = df_h20['RBP_name']
df_h20

fig,ax=plt.subplots(figsize=(10, 15))
sns.heatmap(pd.DataFrame(df_h20['fimo1_motif']), linewidths=0.5, ax=ax, cmap="YlGnBu")
savefn = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/F2.dynamic_region.fimo.enrich.occurance.pdf'
plt.tight_layout()
plt.savefig(savefn)

fig,ax=plt.subplots(figsize=(10, 15))
sns.heatmap(pd.DataFrame(df_h20['odd']), linewidths=0.5, ax=ax, cmap="YlGnBu")
savefn = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/F2.dynamic_region.fimo.enrich.odd.pdf'
plt.tight_layout()
plt.savefig(savefn)

fig,ax=plt.subplots(figsize=(10, 15))
sns.heatmap(pd.DataFrame(df_h20['-log10(pvalue_adj)']), linewidths=0.5, ax=ax, cmap="YlGnBu")
savefn = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/F2.dynamic_region.fimo.enrich.pvalue.pdf'
plt.tight_layout()
plt.savefig(savefn)