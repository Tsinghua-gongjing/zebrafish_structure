import pandas as pd
from nested_dict import nested_dict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
sns.set(style="ticks")
sns.set_context("poster")

def load_upset_file(txt=None, header_ls=None,):
	if txt is None:
		txt = '/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/new_mergepeaks_id_d10/venn.txt'
	set_category_stat_dict = nested_dict(2, int)
	df = pd.read_csv(txt, header=0, sep='\t')
	for sample in df.columns[0:-2]:
		for label,count,name in zip(df[sample], df['Total'], df['Name']):
			if label == 'X' and sample in name:
				set_category_stat_dict[sample][len(name.split('|'))] += count
	set_category_stat_df = pd.DataFrame.from_dict(set_category_stat_dict, orient='index')
	set_category_stat_df = set_category_stat_df.loc[df.columns[0:-2], :]
	print set_category_stat_df

	set_category_stat_df_ratio = set_category_stat_df.div(set_category_stat_df.sum(axis=1), axis=0)
	print set_category_stat_df_ratio

	fig, ax = plt.subplots(1, 2)
	set_category_stat_df.plot(kind='bar', stacked=True, ax=ax[0])
	set_category_stat_df_ratio.plot(kind='bar', stacked=True, ax=ax[1])
	plt.tight_layout()
	savefn = txt.replace('.txt', '.ratio.pdf')
	plt.savefig(savefn)
	plt.close()

def main():
	load_upset_file()

if __name__ == '__main__':
	main()