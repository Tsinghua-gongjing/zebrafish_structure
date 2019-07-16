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
plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/0-figures-format/denovo_motif_scatter_weblogo.sphere_shield.pdf')