import gj,os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="darkgrid")
sns.set_context("poster")
import sys
from nested_dict import nested_dict
import pandas as pd
import numpy as np
from pyfasta import Fasta


def read_fa(fa='/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.fa'):
	gj.printFuncRun('read_fa')
	gj.printFuncArgs()
	fa_dict1 = Fasta(fa, key_fn=lambda key:key.split("\t")[0])
	fa_dict = {i.split()[0]:j[0:] for i,j in fa_dict1.items()}
	print fa_dict.keys()[0:3]
	gj.printFuncRun('read_fa')
	return fa_dict

def read_tmp_out(tmp_out=None,file_str=None,sample=None):
	gj.printFuncRun('read_tmp_out')
	gj.printFuncArgs()
	fa_dict = read_fa()
	tx_base_pos_dict = nested_dict(2, list) # {tx:{'A':[pos1,pos2],'T':[]}}
	base_enrich_dict = nested_dict(1, int)
	with open(tmp_out, 'r') as TMP_OUT:
		for line in TMP_OUT:
			line = line.strip()
			if not line or line.startswith('#'): continue
			arr = line.split('\t')
			transcript_id = arr[0]
			transcript_len = int(arr[1])
			if transcript_len != len(fa_dict[transcript_id]):
				print "transcirpt length not conistent with reference: %s, tmp_out len: %s, reference len: %s"%(transcript_id, transcript_len, len(fa_dict[transcript_id]))
				sys.exit()
			for n,base_enrichment_score in enumerate(arr[4:]):
				score = base_enrichment_score.split(',')[0]
				#if score != "NULL" and float(score) != 0 and float(score) >= 0.3:
				if score != "NULL" and float(score) != 0:
					base = fa_dict[transcript_id][n]
					tx_base_pos_dict[transcript_id][base].append(n)
					base_enrich_dict[base.upper()] += 1
	print base_enrich_dict

	#val_ls = [base_enrich_dict[i] for i in ['A','T','C','G']]
	#gj.plot_ls_pie(labels=['A','T','C','G'],val=val_ls,dic="",title_str="",file_str=file_str)
	TXT = open(file_str, 'w')
	for i,j in base_enrich_dict.items():
		print >>TXT,i+'\t'+str(j)
	TXT.close()

	gj.printFuncRun('read_tmp_out')

def file_info(file_dir=None, result_dir=None):
	if file_dir is None:
		file_dir = '/Share/home/zhangqf7/gongjing/zebrafish/data/RT'
	if result_dir is None:
		result_dir = '/Share/home/zhangqf7/gongjing/zebrafish/result/RTCorrelationPairwise'
	#RT_files = ['DMSO_1cell_rep1', 'DMSO_1cell_rep2', 'NAI_1cell_rep1', 'NAI_1cell_rep2']
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

def main():
	file_info_dict = file_info(file_dir='/Share/home/zhangqf7/gongjing/zebrafish/data/tmpout_by_win_new', result_dir='/Share/home/zhangqf7/gongjing/zebrafish/result/ATCG_content_NewWin/check')
	# tmp_out_ls = [i for i in file_info_dict['paths'] if i.endswith('tmp.out')]
	sample_ls = ['egg', '1cell', '4cell', '64cell', 'sphere', 'shield']
	for tmp_out in tmp_out_ls:
		sample = tmp_out.split('/')[-1].split('.')[0]
		if sample in sample_ls:
			read_tmp_out(tmp_out=tmp_out, file_str=file_info_dict['result_dir'] +'/' + sample+'.icshape.enrich.txt', sample=sample)
	stages_df_ls = [] 
	for j,i in zip(sample_ls, [os.path.join(file_info_dict['result_dir'],"%s.icshape.enrich.txt" % sample) for sample in sample_ls]):
		df = pd.read_csv(i, header=None, index_col=0, names=[j], sep="\t")
		df = df.T
		stages_df_ls.append(df)
	stage_df = pd.concat(stages_df_ls, axis=0)
	stage_df.to_csv(os.path.join(file_info_dict['result_dir'],'total.csv'))

if __name__ == '__main__':
	main()

