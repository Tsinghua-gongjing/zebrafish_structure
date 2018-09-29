import paris_dg2frame
import os
from nested_dict import nested_dict
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
sns.set(style="ticks")
sns.set_context("poster")
import networkx as nx 

def read_dir(dir='/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-5-rep-combine/downsampling_N', to_dgframe=0, get_inter_intra=1, read_nx=1, interRRI_norRNA=1, support_read=3):
	fn_ls = os.listdir(dir)
	# print fn_ls

	fn_stat_dict = nested_dict()
	downsampling_N_draw = dir + '.subnetwork.draw.pdf'
	fig,ax=plt.subplots(10,1)
	for n,fn in enumerate(fn_ls):
		print "process: %s"%(fn)
		dfFile = dir + '/' + fn + '/' + '27-DG'
		frameFile = dfFile + '.txt'
		if to_dgframe:
			paris_dg2frame.DG2Frame(dfFile=dfFile, frameFile=frameFile)
		if get_inter_intra:
			inter, intra = 0, 0
			with open(frameFile, 'r') as TXT:
				for line in TXT:
					line = line.strip()
					if not line or line.startswith('#'):
						continue
					arr = line.split('\t')
					if arr[1] == arr[5]:
						intra += 1
					else:
						inter += 1
			fn_stat_dict[fn]['inter'] = inter
			fn_stat_dict[fn]['intra'] = intra
			fn_stat_dict[fn]['all'] = intra + inter
		if read_nx:
			df = pd.read_csv(frameFile, header=0, sep='\t')
			df['type'] = ['intra' if i == j else 'inter' for i,j in zip(df['lchr'], df['rchr'])]
			df_inter_RRI = df[df['type']=='inter']
			nx_inter_RRI = nx.from_pandas_dataframe(df_inter_RRI, 'lchr', 'rchr')
			fn_stat_dict[fn]['uniq RRI']  = len(nx_inter_RRI.edges())
			if interRRI_norRNA:
				df_inter_RRI = df_inter_RRI[(df_inter_RRI['ltype'].isin(['mRNA', 'lncRNA'])) & (df_inter_RRI['rtype'].isin(['mRNA', 'lncRNA']))]
			df_inter_RRI = df_inter_RRI[df_inter_RRI['support']>=support_read]
			nx_inter_RRI = nx.from_pandas_dataframe(df_inter_RRI, 'lchr', 'rchr')
			nx_inter_RRI_info_dict, G_largest = RRI_network_property2(nx_inter_RRI)
			for i,j in nx_inter_RRI_info_dict.items():
				fn_stat_dict[fn][i] = j
			# fn_stat_dict[fn]['uniq RRI']  = len(nx_inter_RRI.edges())
			if n < 10:
				draw_graph(G_largest, ax=ax[n])
	plt.savefig(downsampling_N_draw)
	savefn = dir + '.stat.txt'
	fn_stat_df = pd.DataFrame.from_dict(fn_stat_dict)
	fn_stat_df = fn_stat_df.T
	fn_stat_df['sampling'] = fn_stat_df.index
	print fn_stat_df.head()

	fn_stat_df.to_csv(savefn, header=True, index=False, sep='\t')

	return fn_stat_df


def RRI_network_property2(nx_inter_RRI=None):
	nx_inter_RRI_info = nx.info(nx_inter_RRI)
	nx_inter_RRI_info_dict = {i.split(':')[0]:i.split(':')[1].strip() for i in nx_inter_RRI_info.split('\n')}
	# nx_inter_RRI_info_dict['average_degree_connectivity'] = nx.average_degree_connectivity(nx_inter_RRI)
	nx_inter_RRI_info_dict['average_clustering'] = nx.average_clustering(nx_inter_RRI)
	g_average_path = []
	"""
	for G in nx.connected_component_subgraphs(nx_inter_RRI):
		try:
			g_average_path.append(nx.average_shortest_path_length(G))
		except:
			pass
	"""
	G_largest = max(nx.connected_component_subgraphs(nx_inter_RRI), key=len)
	# nx_inter_RRI_info_dict['min_average_shortest_path_length'] = np.mean(g_average_path)
	nx_inter_RRI_info_dict['max_average_shortest_path_length'] = nx.average_shortest_path_length(G_largest)
	nx_inter_RRI_info_dict['max_subnetwork_nodes'] = len(G_largest.nodes())
	return nx_inter_RRI_info_dict, G_largest

def readDGFrameFile(filename, interRRI_norRNA=1, support_read=3):
	fn_stat_dict = nested_dict()
	inter, intra = 0, 0
	with open(filename, 'r') as TXT:
		for line in TXT:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			arr = line.split('\t')
			if arr[1] == arr[5]:
				intra += 1
			else:
				inter += 1
	fn_stat_dict['inter'] = inter
	fn_stat_dict['intra'] = intra
	fn_stat_dict['all'] = intra + inter

	df = pd.read_csv(filename, header=0, sep='\t')
	df['type'] = ['intra' if i == j else 'inter' for i,j in zip(df['lchr'], df['rchr'])]
	df_inter_RRI = df[df['type']=='inter']
	nx_inter_RRI = nx.from_pandas_dataframe(df_inter_RRI, 'lchr', 'rchr')
	fn_stat_dict['uniq RRI']  = len(nx_inter_RRI.edges())
	if interRRI_norRNA:
		df_inter_RRI = df_inter_RRI[(df_inter_RRI['ltype'].isin(['mRNA', 'lncRNA'])) & (df_inter_RRI['rtype'].isin(['mRNA', 'lncRNA']))]
	df_inter_RRI = df_inter_RRI[df_inter_RRI['support']>=support_read]
	nx_inter_RRI = nx.from_pandas_dataframe(df_inter_RRI, 'lchr', 'rchr')
	nx_inter_RRI_info_dict, G_largest = RRI_network_property2(nx_inter_RRI)
	for i,j in nx_inter_RRI_info_dict.items():
		fn_stat_dict[i] = j
	# fn_stat_df['sampling'] = ''
	fn_stat_df = pd.DataFrame(fn_stat_dict, index=[0])
	return fn_stat_df

def draw_graph(G, ax):
	nx.draw(G, ax=ax)

def main():
	dir_ls = ['/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-%s-rep-combine/downsampling_N'%(i) for i in [1,4]]
	sample_ls = ['egg', '64cell']
	fn_stat_df_ls = []
	for dir,sample in zip(dir_ls, sample_ls):
		fn_stat_df = read_dir(dir=dir, to_dgframe=0, get_inter_intra=1, support_read=2)
		fn_stat_df['sample'] = sample
		fn_stat_df_ls.append(fn_stat_df)
	fn_stat_df_all = pd.concat(fn_stat_df_ls)
	print fn_stat_df_all.head()
	print fn_stat_df_all.shape

	id_vars_ls = ['sampling', 'sample']
	value_vars_ls = ['all', 'inter', 'intra', 'uniq RRI', 'Average degree', 'average_clustering', 'max_average_shortest_path_length', 'max_subnetwork_nodes']
	fn_stat_df_all_melt = pd.melt(fn_stat_df_all[id_vars_ls+value_vars_ls], id_vars=id_vars_ls, value_vars=value_vars_ls, var_name='Type', value_name='RRI number')
	fn_stat_df_all_melt['RRI number'] = fn_stat_df_all_melt['RRI number'].astype('float')
	print fn_stat_df_all_melt.head()
	print fn_stat_df_all_melt.shape
	fn_stat_df_all_melt.to_csv('/Share/home/zhangqf7/gongjing/zebrafish/result/paris_stat/paris_downsampling.stat.egg_64cell.sp2.txt', header=True, index=True, sep='\t')

	real_file_ls = ['/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-%s-rep-combine/27-DG.txt' % i for i in [1,4]]
	real_stat_df_ls = []
	for dgframefile, sample in zip(real_file_ls, sample_ls):
		real_stat_df = readDGFrameFile(dgframefile)
		real_stat_df['sample'] = sample
		real_stat_df_ls.append(real_stat_df)
	real_stat_df_all = pd.concat(real_stat_df_ls)
	real_stat_df_all.to_csv('/Share/home/zhangqf7/gongjing/zebrafish/result/paris_stat/paris_real.stat.egg_64cell.txt', header=True, index=True, sep="\t")


if __name__ == '__main__':
	main()