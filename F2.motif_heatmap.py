

df_h20 = df.sort_values(by='pvalue_adj', ascending=True).head(100)
df_h20 = df_h20[df_h20['odd']>1]
df_h20 = df_h20[df_h20['RBP_name']!='HuR, ELAVL2/3/4']
df_h20['-log10(pvalue_adj)'] = -np.log10(df_h20['pvalue_adj'])
df_h20 = df_h20[df_h20['-log10(pvalue_adj)']>5]
df_h20.index = df_h20['RBP_name']
df_h20

fig,ax=plt.subplots(figsize=(10, 15))
sns.heatmap(pd.DataFrame(df_h20['fimo1_motif']), linewidths=0.5, ax=ax, cmap="YlGnBu")
savefn = '/Users/gongjing/SeafileSyn/Project/zebrafish/0821-gini-TE/fimo.enrich.denovo.occurance.pdf'
plt.tight_layout()
plt.savefig(savefn)

fig,ax=plt.subplots(figsize=(10, 15))
sns.heatmap(pd.DataFrame(df_h20['odd']), linewidths=0.5, ax=ax, cmap="YlGnBu")
savefn = '/Users/gongjing/SeafileSyn/Project/zebrafish/0821-gini-TE/fimo.enrich.denovo.odd.pdf'
plt.tight_layout()
plt.savefig(savefn)

fig,ax=plt.subplots(figsize=(10, 15))
sns.heatmap(pd.DataFrame(df_h20['-log10(pvalue_adj)']), linewidths=0.5, ax=ax, cmap="YlGnBu")
savefn = '/Users/gongjing/SeafileSyn/Project/zebrafish/0821-gini-TE/fimo.enrich.denovo.pvalue.pdf'
plt.tight_layout()
plt.savefig(savefn)