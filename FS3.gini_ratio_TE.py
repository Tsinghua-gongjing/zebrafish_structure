TE_our = './20190415_TE_corr.txt'
df_TE_our = pd.read_csv(TE_our, header=0, sep='\t')
print df_TE_our.head()

df_TE_our['id'] = df_TE_our['transcript']
df_TE_our = df_TE_our[['id', 'log2(mean(TE(control)))']]
df_TE_our['TE_degree'] = pd.qcut(df_TE_our['log2(mean(TE(control)))'], 4, labels=['q0','q1','q2','q3'])
print df_TE_our.head() , df_TE_our.shape

df_gini_TE = pd.merge(df_gini, df_TE_our, on='id')
print df_gini_TE.head(), df_gini_TE.shape
# print df_gini_TE['TE_degree'].value_counts()

fig, ax = plt.subplots(1,2, figsize=(10, 6))
sns.boxplot(x='TE_degree', y='log2(mean(TE(control)))', data=df_gini_TE, ax=ax[0])
sns.boxplot(x='TE_degree', y='1cell/egg', data=df_gini_TE, ax=ax[1])
plt.tight_layout()
plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/0-figures-format/cell1_egg_gini_ratio_vs_1cell_TE_ourdata.jointplot.pdf')


q0_ls = list(df_gini_TE[df_gini_TE['TE_degree']=='q0']['1cell/egg'])
q3_ls = list(df_gini_TE[df_gini_TE['TE_degree']=='q3']['1cell/egg'])
print stats.mannwhitneyu(q0_ls, q3_ls)
print stats.ttest_ind(q0_ls, q3_ls)