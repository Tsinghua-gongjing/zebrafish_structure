import subprocess
import os
from multiprocessing import Pool
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import gj
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
sns.set(style="ticks")
sns.set_context("poster")
from sklearn.metrics import roc_curve, auc
from nested_dict import nested_dict


def file_info():
	color_stages = sns.color_palette('Set1',n_colors=7, desat=0.8)
	my_pal = {'egg':color_stages[0], '1cell': color_stages[1], '4cell': color_stages[2], '64cell': color_stages[3], '1K': color_stages[4], 'sphere':color_stages[5], 'shield':color_stages[6]}
	file_info_dict = {'const_sample8_v2_pl' : '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/const_sample8_v2.pl',
					  'const_test7_v2_pl'   : '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/const_test7_v2.pl',
					  'random_forest_RBP2_py': '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/random_forest_RBP2.py',
					  'random_forest_RBP3_py': '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/random_forest_RBP3.py',
					  'human_293T_vivo_icshape_out': '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/293T.invivo.out',
					  'human_reference_fasta': '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/gencode.v26.transcripts.ics.std.fa',
					  'zebrafish_reference_fasta': '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/danRer10.refSeq.transcriptome.fa',
					  'HuR_motif_meme': '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/HuR_motif.meme',
					  'HuR_uniq_tx': '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/HuR_uniq_trx.txt',
					  'egg': '/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/egg.icshape.w200.s30.T2.t200.out',
					  '1cell': '/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/1cell.icshape.w200.s30.T2.t200.out',
					  '4cell':'/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/4cell.icshape.w200.s30.T2.t200.out',
					  '64cell':'/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/64cell.icshape.w200.s30.T2.t200.out',
					  '1K':'/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/1K.icshape.w200.s30.T2.t200.out',
					  'sphere':'/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/sphere.icshape.w200.s30.T2.t200.out',
					  'shield':'/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out',
					  'cell_ls': ['egg', '1cell', '4cell', '64cell', '1K', 'sphere', 'shield'],
					  'my_pal':my_pal}
	return file_info_dict

def run_predict(human_pval=0.002, human_null_percent=0.9, zebrafish_pval=0.0001, zebrafish_null_percent=0.6, out_dir='/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR', cell='64cell'):
	"""
	perl const_sample8_v2.pl HuR_motif.meme HuR_uniq_trx.txt HuR gencode.v26.transcripts.ics.std.fa 2e-3 293T.invivo.out 0.9
	perl const_test7_v2.pl HuR_motif.meme HuR danRer10.refSeq.transcriptome.fa 1e-4 64cell.icshape.w200.s30.T2.t200.out 0.6
	python random_forest_RBP2.py HuR_pn_sample_data.txt HuR_pred_motif_ics.txt HuR_zbf_random_forest 0 28 28 
	python random_forest_RBP3.py HuR_pn_sample_data.txt HuR_pred_motif_ics.txt HuR_zbf_random_forest 0 28 28
	"""
	file_info_dict = file_info()
	out_dir = out_dir + '/' + 'human_pval%s_nullpct%s_zebrafish_pval%s_nullpct%s_%s'%(human_pval, human_null_percent, zebrafish_pval, zebrafish_null_percent, cell)
	if os.path.isdir(out_dir):
		print "Exist (%s), will removed"%(out_dir)
		subprocess.call(["rm -rf %s; mkdir %s"%(out_dir, out_dir)], shell=True) 
	else:
		print "Not exist (%s), will creat"%(out_dir)
		subprocess.call(["mkdir %s"%(out_dir)], shell=True) 
	subprocess.call(["cd %s; cp %s ./ ; perl %s HuR_motif.meme %s HuR %s %s %s %s"%(out_dir,
															 file_info_dict['HuR_motif_meme'],
															 file_info_dict['const_sample8_v2_pl'], 
															 file_info_dict['HuR_uniq_tx'],
															 file_info_dict['human_reference_fasta'], 
															 human_pval, 
															 file_info_dict['human_293T_vivo_icshape_out'], 
															 human_null_percent)], shell=True)

	subprocess.call(["cd %s; perl %s HuR_motif.meme HuR %s %s %s %s"%(out_dir,
														  file_info_dict['const_test7_v2_pl'],
												  		  # file_info_dict['HuR_motif_meme'],
												  		  file_info_dict['zebrafish_reference_fasta'],
												  		  zebrafish_pval,
												  		  file_info_dict[cell],
												  		  zebrafish_null_percent)], shell=True)

	subprocess.call(["cd %s; python %s HuR_pn_sample_data.txt HuR_pred_motif_ics.txt HuR_zbf_random_forest 0 28 28; cp random_forest.png random_forest.1.png"%(out_dir, file_info_dict['random_forest_RBP2_py'])], shell=True)
	subprocess.call(["cd %s; python %s HuR_pn_sample_data.txt HuR_pred_motif_ics.txt HuR_zbf_random_forest 0 28 28"%(out_dir, file_info_dict['random_forest_RBP3_py'])], shell=True)

def run_predict_batch(run_predict_pipeline=0, get_auc=1):
	human_pval_ls = [0.00001, 0.0001, 0.0002, 0.001, 0.002, 0.01, 0.02]
	human_null_percent_ls = [0.9, 0.7, 0.5]
	zebrafish_pval_ls = [0.00001, 0.0001, 0.0002, 0.001, 0.002, 0.01, 0.02]
	zebrafish_null_percent_ls = [0.8, 0.6, 0.4]
	cell_ls = ['egg', '1cell', '4cell', '64cell', '1K', 'sphere', 'shield']

	args_tuple_ls = []
	out_dir='/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR'
	for human_pval in human_pval_ls:
		for human_null_percent in human_null_percent_ls:
			for zebrafish_pval in zebrafish_pval_ls:
				for zebrafish_null_percent in zebrafish_null_percent_ls:
					for cell in cell_ls:
						args_tuple_ls.append((human_pval, human_null_percent, zebrafish_pval, zebrafish_null_percent, out_dir, cell))
						# roc100_npz_ls.append('/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/human_pval%s_nullpct%s_zebrafish_pval%s_nullpct%s_%s/HuR_zbf_random_forest_100roc.npz'%(human_pval, human_null_percent, zebrafish_pval, zebrafish_null_percent, cell))
	# print args_tuple_ls

	if run_predict_pipeline:
		pool = Pool(40)
		results = pool.map(run_predict_multiple_wrapper, args_tuple_ls)
		pool.close()

	if get_auc:
		savefn = '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/auc_with_parameters.txt'
		SAVEFN = open(savefn, 'w')
		for (human_pval, human_null_percent, zebrafish_pval, zebrafish_null_percent, out_dir, cell) in args_tuple_ls:
			roc100_npz = '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/human_pval%s_nullpct%s_zebrafish_pval%s_nullpct%s_%s/HuR_zbf_random_forest_100roc.npz'%(human_pval, human_null_percent, zebrafish_pval, zebrafish_null_percent, cell)
			if os.path.isfile(roc100_npz):
				roc_data = np.load(roc100_npz)
				tpr = roc_data['arr_1']
				fpr = roc_data['arr_0']
				roc_auc = auc(fpr, tpr)
				print "exist", roc_auc
				print >>SAVEFN, '\t'.join(map(str, [cell, human_pval, human_null_percent, zebrafish_pval, zebrafish_null_percent, roc_auc]))
			else:
				print "Not exist"
		SAVEFN.close()

def auc_parameter_plot(auc='/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/auc_with_parameters.txt'):
	df = pd.read_csv(auc, header=None, sep='\t')
	df.columns = ['cell', 'h_pval', 'h_null', 'z_pval', 'z_null', 'auc']
	print df.head()

	color_stages = sns.color_palette('Set1',n_colors=7, desat=0.8)
	my_pal = {'egg':color_stages[0], '1cell': color_stages[1], '4cell': color_stages[2], '64cell': color_stages[3], '1K': color_stages[4], 'sphere':color_stages[5], 'shield':color_stages[6]}
	cell_ls = ['egg', '1cell', '4cell', '64cell', '1K', 'sphere', 'shield']

	zebrafish_pval_ls = [0.00001, 0.0001, 0.0002, 0.001, 0.002, 0.01, 0.02]
	zebrafish_null_percent_ls = [0.8, 0.6, 0.4]
	auc_parameter_plot_dir = '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/auc_parameter_plot'
	for zebrafish_pval in zebrafish_pval_ls:
		for zebrafish_null_percent in zebrafish_null_percent_ls:
			df_plot = df[(df['z_pval']==zebrafish_pval) & (df['z_null']==zebrafish_null_percent)]
			if df_plot.shape[0] > 0:
				print "plot"
				print df_plot.head()
				g = sns.FacetGrid(df_plot, col='h_pval', row='h_null', sharex=True, sharey=True, margin_titles=True)
				g.map(sns.barplot, 'cell', 'auc', palette=my_pal, order=cell_ls)
				g.set_xticklabels(rotation=90)
				savefn = auc_parameter_plot_dir + '/' + 'auc_dis.z_pval%s_nullpct%s.png'%(zebrafish_pval, zebrafish_null_percent)
				g.savefig(savefn)
				plt.close()


def run_predict_multiple_wrapper(args):
	return run_predict(*args)

def readIc(ic_file):
	"Read icSHAPE Reactivity File"
	print "read: %s"%(ic_file)
	icDict = dict()
	IN = open(ic_file)
	line = IN.readline()
	while line:
		arr = line.strip().split()
		icDict[arr[0]] = arr[3:]
		line = IN.readline()
	IN.close()
	return icDict

def realdata_stat(realdata_ls=None, positive_cutoff=0.5, extend=0):
	file_info_dict = file_info()
	if realdata_ls is None:
		realdata_ls = ['/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/human_pval0.0002_nullpct0.9_zebrafish_pval0.001_nullpct0.6_%s/HuR_zbf_random_forest_100realdata.txt'%(cell) for cell in file_info_dict['cell_ls']]
	cell_ls = list(set([i.split('/')[-2].split('_')[-1] for i in realdata_ls]))
	cell_icshape_dict = {cell:readIc(file_info_dict[cell]) for cell in cell_ls}
	col6_ls = ['tx', 'length', 'start', 'end', 'pval', 'seq', 'score']
	sample_info_dict = nested_dict()
	df_save_ls = []
	df_save_negative_ls = []
	for realdata in realdata_ls:
		if not os.path.isfile(realdata):
			continue
		print "process: %s"%(realdata)
		cell = realdata.split('/')[-2].split('_')[-1]
		df = pd.read_csv(realdata, header=None, sep=' ')
		print df.head()

		df_positive = df[df[39]>positive_cutoff]
		df_negative = df[df[39]<=positive_cutoff]

		sample_info_dict[cell]['transcript'] = len(set(list(df[0])))
		sample_info_dict[cell]['transcript(postive)'] = len(set(list(df_positive[0])))
		sample_info_dict[cell]['transcript(negative)'] = len(set(list(df_negative[0])))
		sample_info_dict[cell]['regions'] = df.shape[0]
		sample_info_dict[cell]['regions(positive)'] = df_positive.shape[0]
		sample_info_dict[cell]['regions(negative)'] = df_negative.shape[0]

		reactivity_ls = [cell_icshape_dict[cell][tx][max(0, int(start)-1-extend):min(int(end)+extend, int(length))] for tx,start,end,length in zip(df_positive[0],df_positive[2],df_positive[3], df_positive[1])]
		df_save = df_positive.iloc[:,[0,1,2,3,4,5,39]]
		df_save.columns = col6_ls
		df_save['reactivity_ls'] = reactivity_ls
		df_save['reactivity_status'] = [0 if 'NULL' in i else 1 for i in reactivity_ls]
		df_save = df_save[df_save['reactivity_status']==1]
		df_save['gini'] = [gini(i, mode='gini') for i in df_save['reactivity_ls']]
		df_save['mean_reactivity'] = [gini(i, mode='mean_reactivity') for i in df_save['reactivity_ls']]
		df_save['cell'] = cell
		df_save['reactivity_ls'] = [','.join(i) for i in df_save['reactivity_ls']]
		savefn = realdata.replace('.txt', '.ic.txt')
		df_save.to_csv(savefn, header=True, index=False, sep='\t')
		df_save_ls.append(df_save)

		reactivity_ls = [cell_icshape_dict[cell][tx][max(0, int(start)-1-extend):min(int(end)+extend, int(length))] for tx,start,end,length in zip(df_negative[0],df_negative[2],df_negative[3], df_negative[1])]
		df_save_negative = df_negative.iloc[:,[0,1,2,3,4,5,39]]
		df_save_negative.columns = col6_ls
		df_save_negative['reactivity_ls'] = reactivity_ls
		df_save_negative['reactivity_status'] = [0 if 'NULL' in i else 1 for i in reactivity_ls]
		df_save_negative = df_save_negative[df_save_negative['reactivity_status']==1]
		df_save_negative['gini'] = [gini(i, mode='gini') for i in df_save_negative['reactivity_ls']]
		df_save_negative['mean_reactivity'] = [gini(i, mode='mean_reactivity') for i in df_save_negative['reactivity_ls']]
		df_save_negative['cell'] = cell
		df_save_negative['reactivity_ls'] = [','.join(i) for i in df_save_negative['reactivity_ls']]
		df_save_negative_ls.append(df_save_negative)


		sample_info_dict[cell]['regions(positive, stucture)'] = df_save.shape[0]
		sample_info_dict[cell]['transcript(postive, structure)'] = len(set(list(df_save['tx'])))

	sample_info_df = pd.DataFrame.from_dict(sample_info_dict, orient='index')
	print sample_info_df

	savefn = '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/predict_stat.txt'
	cols = ['transcript', 'regions', 'transcript(postive)', 'transcript(negative)', 'regions(positive)', 'regions(negative)', 'transcript(postive, structure)', 'regions(positive, stucture)']
	sample_info_df.loc[file_info_dict['cell_ls'], cols].to_csv(savefn, header=True, index=True, sep='\t')

	df_save_all = pd.concat(df_save_ls)
	savefn = '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/predict_RBP_binding_combine.txt'
	df_save_all.to_csv(savefn, header=True, index=False, sep='\t')

	fig,ax=plt.subplots(1,2, sharey=True)
	sns.boxplot(x='cell', y='gini', data=df_save_all, ax=ax[0], palette=file_info_dict['my_pal'])
	sns.boxplot(x='cell', y='mean_reactivity', data=df_save_all, ax=ax[1], palette=file_info_dict['my_pal'])
	ax[0].set_xticklabels(ax[0].xaxis.get_majorticklabels(), rotation=45)
	ax[1].set_xticklabels(ax[1].xaxis.get_majorticklabels(), rotation=45)
	plt.tight_layout()
	plt.savefig(savefn.replace('.txt', '.png'))
	plt.close()

	df_save_all = pd.concat(df_save_negative_ls)
	savefn = '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/predict_RBP_binding_combine.negative.txt'
	df_save_all.to_csv(savefn, header=True, index=False, sep='\t')

def binding_region_comparison(txt=None, compare_cell_ls=None, compare_2=1, decay_genes=None):
	file_info_dict = file_info()
	if txt is None:
		txt='/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/predict_RBP_binding_combine.txt'
	if compare_cell_ls is None:
		compare_cell_ls = file_info_dict['cell_ls']
	cell_icshape_dict = {cell:readIc(file_info_dict[cell]) for cell in compare_cell_ls}

	df = pd.read_csv(txt, header=0, sep='\t')
	print df.head()

	binding_region_dict = nested_dict()
	compare_cell_binding_region_ls_dict = nested_dict(1, list)
	for tx,start,end,gini_index,mean_reactivity,cell in zip(df['tx'], df['start'], df['end'], df['gini'], df['mean_reactivity'], df['cell']):
		binding_region = ','.join(map(str, [tx, start, end]))
		binding_region_dict[binding_region][cell]['gini'] = gini_index # note: here before set as gini, but we have a function called gini below
		binding_region_dict[binding_region][cell]['mean_reactivity'] = mean_reactivity
		compare_cell_binding_region_ls_dict[cell].append(binding_region)
	print binding_region_dict['NM_199702,2103,2109']
	binding_region_dict = binding_region_dict.to_dict()
	compare_cell_binding_region_ls_dict = compare_cell_binding_region_ls_dict.to_dict()

	gj.venn3plot(mode='string',subsets_ls=[set(compare_cell_binding_region_ls_dict[i]) for i in compare_cell_ls],labels_ls=compare_cell_ls,title_str="overlap",save_fn=txt.replace('.txt', '.overlap.png'),axis=None)

	# common binding regions
	savefn = txt.replace('.txt', '.compare.txt')
	SAVEFN = open(savefn, 'w')
	print >>SAVEFN, '\t'.join(['tx', 'start', 'end'] + ['mean_reactivity', 'type', 'cell'])
	if compare_2:
		common = set(compare_cell_binding_region_ls_dict[compare_cell_ls[0]]) & set(compare_cell_binding_region_ls_dict[compare_cell_ls[1]])
		cell1_specific = set(compare_cell_binding_region_ls_dict[compare_cell_ls[0]]) - set(compare_cell_binding_region_ls_dict[compare_cell_ls[1]])
		cell2_specific = set(compare_cell_binding_region_ls_dict[compare_cell_ls[1]]) - set(compare_cell_binding_region_ls_dict[compare_cell_ls[0]]) 
		for region in common:
			print >>SAVEFN, '\t'.join(map(str, region.split(',')+[ binding_region_dict[region][compare_cell_ls[0]]['mean_reactivity'], 'common', compare_cell_ls[0] ]))
			print >>SAVEFN, '\t'.join(map(str, region.split(',')+[ binding_region_dict[region][compare_cell_ls[1]]['mean_reactivity'], 'common', compare_cell_ls[1] ]))
		for region in cell1_specific:
			print >>SAVEFN, '\t'.join(map(str, region.split(',')+[ binding_region_dict[region][compare_cell_ls[0]]['mean_reactivity'], '%s_specific'%(compare_cell_ls[0]), compare_cell_ls[0] ]))
			tx,start,end = region.split(',')
			if cell_icshape_dict[compare_cell_ls[1]].has_key(tx):
				print >>SAVEFN, '\t'.join(map(str, region.split(',')+[ gini(cell_icshape_dict[compare_cell_ls[1]][tx][(int(start)-1):(int(end))], mode='mean_reactivity'), '%s_specific'%(compare_cell_ls[0]), compare_cell_ls[1] ]))
		for region in cell2_specific:
			tx,start,end = region.split(',')
			if cell_icshape_dict[compare_cell_ls[0]].has_key(tx):
				print >>SAVEFN, '\t'.join(map(str, region.split(',')+[ gini(cell_icshape_dict[compare_cell_ls[0]][tx][(int(start)-1):(int(end))], mode='mean_reactivity'), '%s_specific'%(compare_cell_ls[1]), compare_cell_ls[0] ]))
			print >>SAVEFN, '\t'.join(map(str, region.split(',')+[ binding_region_dict[region][compare_cell_ls[1]]['mean_reactivity'], '%s_specific'%(compare_cell_ls[1]), compare_cell_ls[1] ]))
	else:
		for i,j in binding_region_dict.items():
			if all([1 if j.has_key(c) else 0 for c in compare_cell_ls]):
				mean_reactivity_ls = [j[c]['mean_reactivity'] for c in compare_cell_ls]
				print >>SAVEFN, '\t'.join(map(str, i.split(',')+mean_reactivity_ls))
	SAVEFN.close()

	if compare_2:
		x_order = ['%s_specific'%(compare_cell_ls[0]), 'common', '%s_specific'%(compare_cell_ls[1])]
		df_plot = pd.read_csv(savefn, header=0, sep='\t')
		df_plot = df_plot[df_plot['mean_reactivity']>=0]
		if decay_genes:
			df_plot = df_plot[df_plot['tx'].isin(decay_genes)]
		print df_plot.groupby(['type', 'cell']).count()['tx']
		fig,ax=plt.subplots(figsize=(6,8))
		# gj.venn3plot(mode='string',subsets_ls=[set(compare_cell_binding_region_ls_dict[i]) for i in compare_cell_ls],labels_ls=compare_cell_ls,title_str="overlap",save_fn=txt.replace('.txt', '.overlap.png'),axis=ax[0])
		sns.boxplot(x='type', y='mean_reactivity', hue='cell', data=df_plot, ax=ax, palette=file_info_dict['my_pal'], order=x_order, hue_order=compare_cell_ls)
		ax.legend_.remove()
		plt.tight_layout()
		plt.savefig(savefn.replace('.txt', '.png'))
		plt.close()
	else:
		df_plot = pd.read_csv(savefn, header=0, sep='\t')
		fig,ax=plt.subplots()
		print df_plot.loc[0, compare_cell_ls]
		for i in df_plot.index:
			ax.plot(range(0, len(compare_cell_ls)), df_plot.loc[i, compare_cell_ls], color='grey', alpha=0.3, lw=0.3)
		df_plot_mean = df_plot.loc[:, compare_cell_ls].mean(axis=0)
		ax.plot(range(0, len(compare_cell_ls)), df_plot_mean, color='blue')
		ax.set_title('N=%s'%(df_plot.shape[0]))
		# ax.set_xticklabels(compare_cell_ls, rotation=45)
		plt.xticks(range(0, len(compare_cell_ls)), compare_cell_ls, rotation=45)
		plt.tight_layout()
		plt.savefig(savefn.replace('.txt', '.png'))
		plt.close()

def gini(list_of_values,mode='gini',null_pct=1):
    if len(list_of_values) == 0: return -1
    if list_of_values.count('NULL')/float(len(list_of_values)) > null_pct: return -1
    list_of_values = [i for i in list_of_values if i != 'NULL']
    if len(list_of_values) == 0: return -1
    if type(list_of_values[0]) is str:
        list_of_values = map(float,list_of_values)
    if mode == 'mean_reactivity':
        return np.mean(list_of_values)
    if mode == 'gini':
        if sum(list_of_values) == 0: return 0.67
        sorted_list = sorted(list_of_values)
        height, area = 0, 0
        for value in sorted_list:
            height += value
            area += height - value / 2.
        fair_area = height * len(list_of_values) / 2.
        return (fair_area - area) / fair_area

def read_rpkm2(rpkm=None):
	if rpkm is None:
		rpkm = '/Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelationPairwiseNew/RPKM_combine.merge.txt'
	df = pd.read_csv(rpkm, header=0, sep='\t', index_col=0)
	rpkm_dict = df.to_dict(orient='index')
	#print rpkm_dict.keys()[0:10]
	print rpkm_dict['NM_212847']

	return rpkm_dict, df

def decay_gene_ls(cell1='sphere', cell2='shield', FC_cutoff=2):
	rpkm_dict, df = read_rpkm2()
	decay_genes = []
	for i,j in rpkm_dict.items():
		if j['DMSO_%s'%(cell1)] > 0 and j['DMSO_%s'%(cell2)] > 0:
			fold_change =   np.exp2(j['DMSO_%s'%(cell1)]) / np.exp2(j['DMSO_%s'%(cell2)]) 
			if fold_change >= FC_cutoff:
				decay_genes.append(i)
	print "decay_genes: %s"%(len(decay_genes)), decay_genes[0:5]
	return decay_genes

def main():
	# run_predict()
	# run_predict_batch()
	# auc_parameter_plot()
	# realdata_stat()
	# decay_genes = decay_gene_ls()
	binding_region_comparison(compare_cell_ls=['sphere', 'shield'], decay_genes=None)
	# read_rpkm2()
	

if __name__ == '__main__':
	main()