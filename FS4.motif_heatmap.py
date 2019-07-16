
f = '/Users/gongjing/SeafileSyn/Project/zebrafish/0821-gini-TE/enrich_sphere_shield/fimo.enrich.txt'
df = pd.read_csv(f, header=0, sep='\t')
print df.head()

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
plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/0821-gini-TE/enrich_sphere_shield/fimo.enrich.txt.pdf')

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

savefn = '/Users/gongjing/SeafileSyn/Project/zebrafish/0821-gini-TE/enrich_sphere_shield/fimo.enrich.denovo.occurance.pdf'
plt.tight_layout()
plt.savefig(savefn)