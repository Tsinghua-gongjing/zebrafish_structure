gj_gini = '/Users/soul/BaiduyunDisk/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/global_str_change/data/RPKM_combine.merge.addRP33.t200.gini.null40.txt'
df_gj_gini = pd.read_csv(gj_gini,header=0, index_col=0, sep='\t', na_values='NULL')
sample_ls = ['egg', '1cell', '1cell-RK-33', '4cell', '64cell', 'sphere', 'shield']
feature = 'transcript'
col_selected = ['%s(%s)'%(j,feature) for j in sample_ls]
df_gj_gini_tx = df_gj_gini[col_selected]
df_gj_gini_tx.columns = sample_ls

plot_data = df_gj_gini_tx[['egg', '1cell']].dropna(how='all')
plot_data['id'] = plot_data.index
# plot_data = pd.melt(plot_data, id_vars=['id'], value_vars=['egg', '1cell'], value_name='Gini index')
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(7,6))
sns.distplot(plot_data['egg'].dropna(), bins=100,kde=False, norm_hist=0, color='#d02e30',ax=ax, label='Egg')
sns.distplot(plot_data['1cell'].dropna(), bins=100,kde=False, norm_hist=0, color='#447dab', ax=ax, label='1 Cell')
ax.axvline(np.median(plot_data['egg'].dropna()), linestyle="--", color='#d02e30')
ax.axvline(np.median(plot_data['1cell'].dropna()), linestyle="--", color='#447dab')
sns.despine(top=True, right=True)
plt.xlabel('Gini index')
plt.ylabel('Number of transcript')
# plt.xlim((0.3,0.7))
plt.legend()
# plt.savefig('../icSHAPE/1cell_change/comparison_gini_1cell_vs_egg.pdf', dpi=220)
print 'egg', np.median(plot_data['egg'].dropna())
print '1cell', np.median(plot_data['1cell'].dropna())
# stats.ttest_ind(plot_data['egg'].dropna(), plot_data['1cell'].dropna(), )
stats.mannwhitneyu(plot_data['egg'].dropna(), plot_data['1cell'].dropna(), alternative='greater')

plot_data = df_gj_gini_tx[['1cell', '1cell-RK-33']].dropna(how='all')
plot_data['id'] = plot_data.index
# plot_data = pd.melt(plot_data, id_vars=['id'], value_vars=['egg', '1cell'], value_name='Gini index')
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(7,6))
sns.distplot(plot_data['1cell-RK-33'].dropna(), bins=100,kde=False, norm_hist=0, color='#d02e30',ax=ax, label='1 Cell(RK-33)')
sns.distplot(plot_data['1cell'].dropna(), bins=100,kde=False, norm_hist=0, color='#447dab', ax=ax, label='1 Cell')
ax.axvline(np.median(plot_data['1cell-RK-33'].dropna()), linestyle="--", color='#d02e30')
ax.axvline(np.median(plot_data['1cell'].dropna()), linestyle="--", color='#447dab')
sns.despine(top=True, right=True)
plt.xlabel('Gini index')
plt.ylabel('Number of transcript')
# plt.xlim((0.3,0.7))
plt.legend()
plt.savefig('../icSHAPE/1cell_change/cell1_RK33_gini_hist.pdf', dpi=220)
print '1cell-RK-33', np.median(plot_data['1cell-RK-33'].dropna())
print '1cell', np.median(plot_data['1cell'].dropna())

# stats.ttest_ind(plot_data['1cell-RK-33'].dropna(), plot_data['1cell'].dropna())
stats.mannwhitneyu(plot_data['1cell-RK-33'].dropna(), plot_data['1cell'].dropna(), alternative='greater')