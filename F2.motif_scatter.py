from adjustText import adjust_text

f = '/Users/gongjing/SeafileSyn/Project/zebrafish/0821-gini-TE/fimo.enrich.txt'
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
plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/0821-gini-TE/fimo.enrich.txt.pdf')