# TE_polyA.ipynb

df_gini = pd.read_csv('./RPKM_combine.merge.addRP33.t200.gini.null40.txt', header=0, index_col=0, sep='\t')
df_gini = df_gini[['1cell-RK-33(transcript)', '1cell(transcript)', 'egg(transcript)']]
# df_gini = df_gini[['1cell-RK-33(transcript)']]
df_gini.dropna(inplace=True)
# df_gini['open_degree'] = pd.qcut(df_gini['1cell-RK-33(transcript)'], 4, labels=['q0', 'q1', 'q2', 'q3'])
print df_gini.head(), df_gini.shape

df_gini['RK33/1cell'] = df_gini['1cell-RK-33(transcript)'] / df_gini['1cell(transcript)']
df_gini['egg/1cell'] = df_gini['egg(transcript)'] / df_gini['1cell(transcript)']
df_gini.head()

sns.jointplot(x='RK33/1cell', y='egg/1cell', data=df_gini, kind="reg", stat_func=stats.pearsonr, size=10)
plt.savefig('./egg_1cell_RK33_corr.pdf')