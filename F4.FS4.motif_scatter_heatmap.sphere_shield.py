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



site_num_ls2 = [16132,9885,4869,2019,1974,251,924,168,534,374,298,485,545]
pvalue_ls2 = [1.8e-100,2.6e-80,1.4e-066,2.8e-028,5.6e-023,1.6e-016,2.4e-015,8.1e-011,2.9e-010,4.7e-009,1.8e-008,1.3e-007,7.6e-007]

dd_sphere_shield = pd.DataFrame({'site_num':site_num_ls2, 'pvalue':pvalue_ls2})
dd_sphere_shield['-log10(pvalue)'] = -np.log10(dd_sphere_shield['pvalue'])
dd_sphere_shield['log2(site_num)'] = np.log2(dd_sphere_shield['site_num'])
dd_sphere_shield

fig,ax = plt.subplots(figsize=(10,10))
dd_sphere_shield.plot(kind='scatter', x='log2(site_num)', y='-log10(pvalue)', ax=ax)
ax.set_xlim(0,)
plt.tight_layout()
plt.savefig('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/F4.dynamic_region.sphere_shield.UTR3.fimo.enrich.pdf')





f = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/dynamic_region/sphere_shield/fimo.enrich.txt'
df = pd.read_csv(f, header=0, sep='\t')
print df.head()

df['pvalue_adj'] = [10e-300 if float(i)==0 else i for i in df['pvalue_adj']]
df = df[df['pvalue_adj']<=0.05]
df['log2(fimo1_motif)'] = np.log2(df['fimo1_motif'])
df['-log10(pvalue_adj)'] = -np.log10(df['pvalue_adj'])

# print df

# fig,ax = plt.subplots(figsize=(25,20))
# df.plot(kind='scatter', x='log2(fimo1_motif)', y='-log10(pvalue_adj)', ax=ax)
# ax.set_xlim(0,)

# texts = []
# for x,y,t in zip(df['log2(fimo1_motif)'], df['-log10(pvalue_adj)'], df['RBP_name']):
# #     print x,y,t
#     texts.append(plt.text(x, y, t, fontsize=20))
# adjust_text(texts, only_move={'text': 'xy'})

# plt.tight_layout()
# plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/0821-gini-TE/enrich_sphere_shield/fimo.enrich.txt.pdf')

df_h20 = df.sort_values(by='pvalue_adj', ascending=True).head(100)
df_h20 = df_h20[df_h20['odd']>1]
df_h20 = df_h20[df_h20['RBP_name']!='HuR, ELAVL2/3/4']
df_h20['-log10(pvalue_adj)'] = -np.log10(df_h20['pvalue_adj'])
df_h20 = df_h20[df_h20['-log10(pvalue_adj)']>5]
df_h20.index = df_h20['RBP_name']
df_h20

fig,ax=plt.subplots(1,3, figsize=(30, 15))
sns.heatmap(pd.DataFrame(df_h20['fimo1_motif']), linewidths=0.5, cmap="YlGnBu", ax=ax[0])
sns.heatmap(pd.DataFrame(df_h20['odd']), linewidths=0.5, cmap="YlGnBu", ax=ax[1])
sns.heatmap(pd.DataFrame(df_h20['-log10(pvalue_adj)']), linewidths=0.5, cmap="YlGnBu", ax=ax[2])

savefn = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/FS4.fimo.enrich.heatmap.pdf'
plt.tight_layout()
plt.savefig(savefn)