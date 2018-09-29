import pandas as pd
import numpy as np
import matplotlib 
import seaborn as sns
from matplotlib import pyplot as plt
import os

netstats_file = '/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/PARIS/network_stats/paris_downsampling.stat.egg_64cell.sp2.txt'
fn_stat_df_all_melt = pd.read_csv(netstats_file, header=0, index_col=0, sep="\t")
print fn_stat_df_all_melt.head()
real_stat_file = '/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/PARIS/network_stats/paris_real.stat.egg_64cell.txt'
real_stat_df_all = pd.read_csv(real_stat_file, header=0, index_col=0, sep="\t")
real_stat_df_all.index = real_stat_df_all['sample']
print real_stat_df_all.head()

color_stages = sns.color_palette('Set1',n_colors=7, desat=0.8)
my_pal = {'egg':color_stages[0], '64cell': color_stages[3], }
sns.set_style('ticks')

fig, ax = plt.subplots(figsize=(6,6))
fn_stat_df_all_melt_RRI_property = fn_stat_df_all_melt[fn_stat_df_all_melt['Type'].isin(['Average degree', 'average_clustering', 'max_average_shortest_path_length', 'max_subnetwork_nodes'])]
sns.boxplot(y='RRI number', x='sample', data=fn_stat_df_all_melt_RRI_property[fn_stat_df_all_melt_RRI_property['Type']=='average_clustering'], ax=ax, palette=my_pal)
ax.set_ylabel('Average clustering coefficient')
sample_ls = ['egg', '64cell']
color_map = [my_pal[sam] for sam in sample_ls]
pos = range(len(sample_ls))
ax.scatter(y=real_stat_df_all['average_clustering'], x=pos, c=color_map, s=120, edgecolors='white', linewidth=2)
plt.tick_params(top=False, right=False)
plt.xlabel('')
plt.savefig('average_clustering.egg_64cell.downsampling.pdf', dpi=200, bbox_inches='tight')

sns.set_style('ticks')
color_stages = sns.color_palette('Set1',n_colors=7, desat=0.8)
my_pal = {'egg':color_stages[0], '4cell': color_stages[2], '64cell': color_stages[3], '1K': color_stages[4]}
fig, ax = plt.subplots(figsize=(6,6))
fn_stat_df_all_melt_RRI_property = fn_stat_df_all_melt[fn_stat_df_all_melt['Type'].isin(['Average degree', 'average_clustering', 'max_average_shortest_path_length', 'max_subnetwork_nodes'])]
sns.boxplot(y='RRI number', x='sample', data=fn_stat_df_all_melt_RRI_property[fn_stat_df_all_melt_RRI_property['Type']=='Average degree'], ax=ax, palette=my_pal)
ax.set_ylabel('Average degree')
sample_ls = ['egg', '64cell']
color_map = [my_pal[sam] for sam in sample_ls]
pos = range(len(sample_ls))
ax.scatter(y=real_stat_df_all['Average degree'], x=pos, c=color_map, s=120, edgecolors='white', linewidth=2)
plt.tick_params(top=False, right=False)
plt.xlabel('')
plt.savefig('average_degree.egg_64cell.downsampling.pdf', dpi=200, bbox_inches='tight')