# /Share2/home/zhangqf7/gongjing/zebrafish/data/RPKM_new
# [zhangqf7@ZIO01 RPKM_new]$ grep NM_131452 DMSO*
# DMSO_1cell_rep1:NM_131452       2133    73155   0       269.213107543935
# DMSO_1cell_rep2:NM_131452       2133    43305   0       270.746439888154
# DMSO_1K_rep1:NM_131452  2133    75425   0       362.014081473829
# DMSO_1K_rep2:NM_131452  2133    68142   0       232.019760230161
# DMSO_4cell_rep1:NM_131452       2133    50131   0       261.870572356349
# DMSO_4cell_rep2:NM_131452       2133    55652   0       279.432526454006
# DMSO_64cell_rep1:NM_131452      2133    46609   0       200.216480594238
# DMSO_64cell_rep2:NM_131452      2133    45895   0       184.101487792659
# DMSO_egg_rep1:NM_131452 2133    52011   0       268.665425072292
# DMSO_egg_rep2:NM_131452 2133    69552   0       261.202857921502
# DMSO_epiboly_rep1:NM_131452     2133    60974   0       393.294976408747
# DMSO_epiboly_rep2:NM_131452     2133    101602  0       369.526269255225
# DMSO_shield_rep1:NM_131452      2133    97726   0       479.140603155748
# DMSO_shield_rep2:NM_131452      2133    128909  0       500.724656957156
# DMSO_sphere_rep1:NM_131452      2133    82802   0       277.643785739012
# DMSO_sphere_rep2:NM_131452      2133    82447   0       256.748858470424

# [zhangqf7@ZIO01 RPKM_new]$ grep NM_130909 DMSO*
# DMSO_1cell_rep1:NM_130909       2271    2843    0       9.82658844538043
# DMSO_1cell_rep2:NM_130909       2271    1685    0       9.89460218994131
# DMSO_1K_rep1:NM_130909  2271    683     0       3.07896370845702
# DMSO_1K_rep2:NM_130909  2271    887     0       2.83666100492287
# DMSO_4cell_rep1:NM_130909       2271    1545    0       7.58023264069506
# DMSO_4cell_rep2:NM_130909       2271    2093    0       9.87049949208604
# DMSO_64cell_rep1:NM_130909      2271    788     0       3.17928877416371
# DMSO_64cell_rep2:NM_130909      2271    1495    0       5.63257311515208
# DMSO_egg_rep1:NM_130909 2271    1729    0       8.38851896369865
# DMSO_egg_rep2:NM_130909 2271    3235    0       11.4108056320287
# DMSO_epiboly_rep1:NM_130909     2271    4121    0       24.9660622183875
# DMSO_epiboly_rep2:NM_130909     2271    6478    0       22.1287929517002
# DMSO_shield_rep1:NM_130909      2271    3801    0       17.5034816235738
# DMSO_shield_rep2:NM_130909      2271    4749    0       17.3257334672649
# DMSO_sphere_rep1:NM_130909      2271    946     0       2.97928457606504
# DMSO_sphere_rep2:NM_130909      2271    892     0       2.60898894798169

elavl1a = pd.DataFrame({'egg':[268.665425072292,261.202857921502], 
                  '1cell':[269.213107543935,270.746439888154],
                  '4cell':[261.870572356349,279.432526454006],
                  '64cell':[200.216480594238,184.101487792659],
                  'sphere':[277.643785739012,256.748858470424],
                  'shield':[479.140603155748,500.724656957156]})


elavl1b = pd.DataFrame({'egg':[8.38851896369865,11.4108056320287], 
                  '1cell':[9.82658844538043,9.89460218994131],
                  '4cell':[7.58023264069506,9.87049949208604],
                  '64cell':[3.17928877416371,5.63257311515208],
                  'sphere':[2.97928457606504,2.60898894798169],
                  'shield':[17.5034816235738,17.3257334672649]})

sample_ls = ['egg', '1cell', '4cell', '64cell', 'sphere', 'shield']
elavl1a = elavl1a.loc[:, sample_ls]
elavl1b = elavl1b.loc[:, sample_ls]

elavl1a_mean = elavl1a.mean()
elavl1a_std = elavl1a.std()
print elavl1a_mean
print elavl1a_std

elavl1b_mean = elavl1b.mean()
elavl1b_std = elavl1b.std()
print elavl1b_mean
print elavl1b_std

ind = np.arange(len(elavl1a_mean))
width = 0.35


fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, elavl1a_mean, width, yerr=elavl1a_std,
                label='Elavl1a')
rects2 = ax.bar(ind + width/2, elavl1b_mean, width, yerr=elavl1b_std,
                label='Elavl1b')
ax.legend()

plt.tight_layout()
plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/RPKM_elavl1.errbar.pdf')