import subprocess
import itertools
from nested_dict import nested_dict
import sys,os
import pandas as pd, numpy as np

def file_info(file_dir=None, result_dir=None):
	if file_dir is None:
		file_dir = '/Share/home/zhangqf7/gongjing/zebrafish/data/RT'
	if result_dir is None:
		result_dir = '/Share/home/zhangqf7/gongjing/zebrafish/result/RTCorrelationPairwise'
	files = os.listdir(file_dir)
	NAI_files = [i for i in files if i.startswith('NAI')]
	DMSO_files = [i for i in files if i.startswith('DMSO')]
	paths = [file_dir+'/'+i for i in files]
	file_info_dict = nested_dict()

	file_info_dict['file_dir'] = file_dir
	file_info_dict['files'] = files
	file_info_dict['paths'] = paths
	file_info_dict['result_dir'] = result_dir

	return file_info_dict.to_dict()

def read_rpkm(rpkm_txt=None):
	if rpkm_txt is None:
		rpkm_txt = '/Share/home/zhangqf7/gongjing/zebrafish/data/RPKM/DMSO_1cell_rep1'
	rpkm_dict = nested_dict()
	with open(rpkm_txt, 'r') as TXT:
		for line in TXT:
			line = line.strip()
			if not line: continue
			if line.startswith('#'):
				header = line.replace('#', '').split('\t')
				continue
			arr = line.split('\t')
			for i,j in zip(header, arr):
				rpkm_dict[arr[0]][i] = j 
	return rpkm_dict.to_dict()

def compare_correlation_rpkm(savefn=None, file_dir='/Share/home/zhangqf7/gongjing/zebrafish/data/RPKM', result_dir='/Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelationPairwise'):
	if savefn is None:
		savefn = result_dir + '/'+ 'RPKM_combine.txt' #'/Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelation/RPKM_combine.txt'
	
	rpkm_file_info_dict = file_info(file_dir=file_dir, result_dir=result_dir)
	rpkm_file_rpkmdict_dict = {}
	sample_ls = []
	for rpkm_txt in rpkm_file_info_dict['paths']:
		sample = rpkm_txt.split('/')[-1]
		rpkm_file_rpkmdict_dict[sample] = read_rpkm(rpkm_txt=rpkm_txt)
		sample_ls.append(sample)

	sample_compare_ls = itertools.combinations_with_replacement(sample_ls,2)
	rpkm_pearson_df_ls = []
	for (sample1, sample2) in sample_compare_ls:
		print "compare: %s, %s"%(sample1, sample2)

		savefn = rpkm_file_info_dict['result_dir'] + '/' + sample1+'-'+sample2
		rpkm_pearson = compare_correlation_rpkm_pairwise(rpkm_dict1=rpkm_file_rpkmdict_dict[sample1], rpkm_dict2=rpkm_file_rpkmdict_dict[sample2], savefn=savefn, sample_ls=[sample1, sample2])

		rpkm_pearson_df = pd.DataFrame({'sample1':[sample1], 'sample2':[sample2], 'pearson':[rpkm_pearson]})
		print rpkm_pearson_df

		rpkm_pearson_df_ls.append(rpkm_pearson_df)
	savefn = rpkm_file_info_dict['result_dir'] + '/' + 'stat.spearman.txt'
	rpkm_pearson_df_ls_combine = pd.concat(rpkm_pearson_df_ls)
	rpkm_pearson_df_ls_combine.to_csv(savefn, header=True, index=False, sep='\t')

	sample_tx_ls = [j.keys() for i,j in rpkm_file_rpkmdict_dict.items()]
	sample_tx_overlap = list(set(sample_tx_ls[0]).intersection(*sample_tx_ls))
	print len(sample_tx_ls), [len(i) for i in sample_tx_ls], len(sample_tx_overlap)

	sample_ls = [i.split('/')[-1] for i in rpkm_file_info_dict['paths']]
	print sample_ls

	with open(result_dir + '/'+ 'RPKM_combine.txt', 'w') as TXT:
		print >>TXT, '\t'.join(sample_ls)
		for tx in sample_tx_overlap:
			tx_rpkm_str = '\t'.join([rpkm_file_rpkmdict_dict[i][tx]['RPKM'] for i in sample_ls])
			print >>TXT, tx+'\t'+tx_rpkm_str

def compare_correlation_rpkm_pairwise(rpkm_dict1=None, rpkm_dict2=None, savefn=None, sample_ls=None, min_rpkm=1):
	sample_tx_overlap = list(set(rpkm_dict1.keys()) & set(rpkm_dict2.keys()))
	with open(savefn, 'w') as TXT:
		print >>TXT, '\t'.join(sample_ls)+'\t'+'transcript'
		for tx in sample_tx_overlap:
			if float(rpkm_dict1[tx]['RPKM']) >= min_rpkm and float(rpkm_dict2[tx]['RPKM']) >= min_rpkm:
				print >>TXT, rpkm_dict1[tx]['RPKM'] +'\t' +rpkm_dict2[tx]['RPKM'] + '\t' + tx
	df = pd.read_csv(savefn, header=0, sep='\t')
	df_similarity_pearson = df[sample_ls].corr(method='spearman')

	return df_similarity_pearson.iloc[0,1]

def main():
	compare_correlation_rpkm(file_dir='/Share/home/zhangqf7/gongjing/zebrafish/data/RPKM_new', result_dir='/Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelationPairwiseNew')

if __name__ == '__main__':
	main()
