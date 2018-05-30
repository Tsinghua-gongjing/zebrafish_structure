import subprocess
from nested_dict import nested_dict
import pandas as pd 
import numpy as np

def file_info():
	"""
	load related file info into dict
	"""
	file_info_dict = nested_dict()

	data_dir = '/Share/home/zhangqf7/gongjing/zebrafish/data/paris'
	sample_num_ls = ['shi-zp-1-rep-combine', 'shi-zp-2-rep-combine', 'shi-zp-3-rep-combine', 'shi-zp-4-rep-combine',
					 'shi-zp-5-rep-combine', 'shi-zp-7-rep-combine', 'shi-zp-8-rep-combine']
	sample_ls = ['egg', '1cell', '4cell', '64cell', '1k', 'shield', 'epiboly']

	file_info_dict['sample_num_ls'] = sample_num_ls
	file_info_dict['sample_ls'] = sample_ls

	for sample_num, sample in zip(sample_num_ls, sample_ls):
		sample_dir = data_dir + '/' + sample_num
		DG = sample_dir + '/' + '27-DG'
		DG_txt = sample_dir + '/' + '27-DG.txt'
		star_log = sample_dir + '/' + '6-starLog.final.out'
		fq_ori = sample_dir + '/' + 'clean_data/1-ori.fastq'
		fq_cut = sample_dir + '/' + 'clean_data/2-cut.fastq'
		fq_cut_filter = sample_dir + '/' + 'clean_data/2-cut_filter.fastq'
		fq_derep = sample_dir + '/' + 'clean_data/3-derep.fastq'
		fq_trim = sample_dir + '/' + 'clean_data/4-trim.fastq'
		fq_qua = sample_dir + '/' + 'clean_data/5-qua.fastq'
		align_out_sam = sample_dir + '/' + '6-starAligned.out.sam'
		chimeric_out_sam = sample_dir + '/' + '6-starChimeric.out.sam'

		file_info_dict[sample]['sample_dir'] = sample_dir
		file_info_dict[sample]['DG'] = DG
		file_info_dict[sample]['DG_txt'] = DG_txt
		file_info_dict[sample]['star_log'] = star_log
		file_info_dict[sample]['fq_ori'] = fq_ori
		file_info_dict[sample]['fq_cut'] = fq_cut
		file_info_dict[sample]['fq_cut_filter'] = fq_cut_filter
		file_info_dict[sample]['fq_derep'] = fq_derep
		file_info_dict[sample]['fq_trim'] = fq_trim
		file_info_dict[sample]['fq_qua'] = fq_qua
		file_info_dict[sample]['align_out_sam'] = align_out_sam
		file_info_dict[sample]['chimeric_out_sam']= chimeric_out_sam

	return file_info_dict.to_dict()

def sample_RRI_type_dis():
	"""
	output statistics of inter RRI type distribution

	[zhangqf7@loginview02 shi-zp-1-rep-combine]$ cat 27-DG.inter.stat.txt
	data 1 2 3 4 5 6 7 8 9
	data 202,75,78 83,169,102 205,185,111 98,180,208 129,112,182 238,130,238 255,140,0 74,113,178 169,169,169
	data mRNA lncRNA miRNA misc_RNA pseudogene rRNA snRNA snoRNA other
	mRNA 1017.0 0.0 nan 0.0 0.0 0.0 0.0 0.0 0.0
	lncRNA 201.0 15.0 nan 0.0 0.0 67.0 0.0 0.0 1.0
	miRNA 6.0 0.0 nan 0.0 0.0 0.0 0.0 0.0 0.0
	misc_RNA 19.0 4.0 nan 0.0 0.0 0.0 0.0 0.0 0.0
	pseudogene 7.0 0.0 nan 0.0 0.0 1.0 0.0 0.0 0.0
	rRNA 3403.0 468.0 nan 139.0 2.0 1302.0 8.0 10.0 74.0
	snRNA 11.0 3.0 nan 1.0 0.0 26.0 0.0 0.0 1.0
	snoRNA 4.0 0.0 nan 0.0 0.0 4.0 0.0 1.0 0.0
	other 68.0 4.0 nan 0.0 0.0 55.0 0.0 0.0 1.0


	"""
	file_info_dict = file_info()
	sample_ls = ['egg', '1cell', '4cell', '64cell', '1k']
	RNA_type_ls = ['mRNA', 'lncRNA', 'miRNA', 'misc_RNA', 'pseudogene', 'rRNA', 'snRNA', 'snoRNA', 'other']
	RNA_type_color_ls = ['202,75,78', '83,169,102', '205,185,111', '98,180,208', '129,112,182', '238,130,238', '255,140,0', '74,113,178', '169,169,169']
	RNA_type_color_dict = {i:j for i,j in zip(RNA_type_ls, RNA_type_color_ls)}
	sample_inter_intra_dict = nested_dict()
	for sample in sample_ls:
		DG_txt = file_info_dict[sample]['DG_txt']
		df = pd.read_csv(DG_txt, header=0, sep='\t')
		df = df[df['support'] >= 3]

		RRI_type = ['intra' if i== j else 'inter' for i,j in zip(df['lchr'], df['rchr'])]
		df['RRI_type'] = RRI_type
		#print df.head()

		df_inter = df[df['RRI_type']=='inter']
		df_intra = df[df['RRI_type']=='intra']

		df_inter_stat = df_inter.groupby(['ltype', 'rtype']).count().unstack().fillna(0)['#Group']
		df_inter_stat = df_inter_stat.loc[RNA_type_ls, RNA_type_ls]
		print df_inter_stat

		print 'iner', df_inter.shape
		print 'intra', df_intra.shape
		sample_inter_intra_dict[sample]['inter-RRI'] = df_inter.shape[0]
		sample_inter_intra_dict[sample]['inter-RNA'] = len(set(list(df_inter['lchr']) + list(df_inter['rchr'])))
		sample_inter_intra_dict[sample]['intra-RRI'] = df_intra.shape[0]
		sample_inter_intra_dict[sample]['intra-RNA'] = len(set(list(df_intra['lchr']) + list(df_intra['rchr'])))
		sample_inter_intra_dict[sample]['RRI'] = df_inter.shape[0] + df_intra.shape[0]

		savefn_inter = DG_txt.replace('.txt', '.inter.stat.txt')
		with open(savefn_inter, 'w') as SAVEFN:
			print >>SAVEFN, ' '.join(map(str, ['data']+range(1, len(RNA_type_ls)+1 )))
			print >>SAVEFN, ' '.join(['data']+RNA_type_color_ls)
			print >>SAVEFN, ' '.join(['data']+RNA_type_ls)
			for idx,row in df_inter_stat.iterrows():
				print >>SAVEFN, ' '.join(map(str, [idx]+list(row)))

	print pd.DataFrame.from_dict(sample_inter_intra_dict, orient='index')


def RRI_type_dist_conversion(txt='/Share/home/zhangqf5/gongjing/DNA-RNA-Protein-interaction-correlation-12-18/data/type_dis/RRI_union_deduplicates.split_full.type_dis.human.txt'):
	"""
	format conversion of type distribution

	from:
	data 1 2 3 4 5 6 7 8 9
	data 202,75,78 83,169,102 205,185,111 98,180,208 129,112,182 238,130,238 255,140,0 74,113,178 169,169,169
	data mRNA lncRNA miRNA misc_RNA pseudogene rRNA snRNA snoRNA other
	mRNA 1017.0 0.0 nan 0.0 0.0 0.0 0.0 0.0 0.0
	lncRNA 201.0 15.0 nan 0.0 0.0 67.0 0.0 0.0 1.0
	miRNA 6.0 0.0 nan 0.0 0.0 0.0 0.0 0.0 0.0
	misc_RNA 19.0 4.0 nan 0.0 0.0 0.0 0.0 0.0 0.0
	pseudogene 7.0 0.0 nan 0.0 0.0 1.0 0.0 0.0 0.0
	rRNA 3403.0 468.0 nan 139.0 2.0 1302.0 8.0 10.0 74.0
	snRNA 11.0 3.0 nan 1.0 0.0 26.0 0.0 0.0 1.0
	snoRNA 4.0 0.0 nan 0.0 0.0 4.0 0.0 1.0 0.0
	other 68.0 4.0 nan 0.0 0.0 55.0 0.0 0.0 1.0

	to:
	data 1 2 3 4 5 6 7 8 9
	data 202,75,78 83,169,102 205,185,111 98,180,208 129,112,182 238,130,238 255,140,0 74,113,178 169,169,169
	data mRNA lncRNA miRNA misc_RNA pseudogene rRNA snRNA snoRNA other
	mRNA 1017.0 201.0 0 19.0 7.0 3403.0 11.0 4.0 68.0
	lncRNA 0 15.0 0 4.0 0.0 535.0 3.0 0.0 5.0
	miRNA 0 0 0 0 0 0 0 0 0
	misc_RNA 0 0 0 0.0 0.0 139.0 1.0 0.0 0.0
	pseudogene 0 0 0 0 0.0 3.0 0.0 0.0 0.0
	rRNA 0 0 0 0 0 1302.0 34.0 14.0 129.0
	snRNA 0 0 0 0 0 0 0.0 0.0 1.0
	snoRNA 0 0 0 0 0 0 0 1.0 0.0
	other 0 0 0 0 0 0 0 0 1.0
	"""
	dis_dict = nested_dict(2, int)
	dis_revise_dict = nested_dict(2, int)
	with open(txt, 'r') as TXT:
		for n,line in enumerate(TXT):
			line = line.strip()
			print n,line
			if n == 0:
				dis_dict['row_order'] = line
			elif n == 1:
				dis_dict['row_color'] = line
			elif n == 2:
				dis_dict['row_rnas'] = line
				col_rna_ls = line.split()[1:]
			else:
				row_rna = line.split()[0]
				val_ls = map(float, line.split()[1:])
				for col_rna, val in zip(col_rna_ls, val_ls):
					dis_dict[row_rna][col_rna] = val
	print dis_dict

	rna_pair_ls = []
	for row_rna in col_rna_ls:
		for col_rna in col_rna_ls:
			dis_revise_dict[row_rna][col_rna] = 0

	for row_rna in col_rna_ls:
		for col_rna in col_rna_ls:
			rna_pair1 = row_rna + '-' + col_rna
			rna_pair2 = col_rna + '-' + row_rna
			if rna_pair1 in rna_pair_ls: 
				continue
			if rna_pair2 in rna_pair_ls:
				continue
			rna_pair_ls.append(rna_pair1)
			rna_pair_ls.append(rna_pair2)
			dis_revise_dict[row_rna][col_rna] += dis_dict[row_rna][col_rna]
			if row_rna == col_rna:
				continue
			dis_revise_dict[row_rna][col_rna] += dis_dict[col_rna][row_rna]
	print dis_revise_dict

	savefn = txt.replace('txt', 'revise.txt')
	with open(savefn, 'w') as SAVEFN:
		print >>SAVEFN, dis_dict['row_order']
		print >>SAVEFN, dis_dict['row_color']
		print >>SAVEFN, dis_dict['row_rnas']
		for row_rna in col_rna_ls:
			row_rna_ls = [dis_revise_dict[row_rna][col_rna] for col_rna in col_rna_ls]
			row_rna_ls = [0 if np.isnan(i) else i for i in row_rna_ls]
			print >>SAVEFN,row_rna+' '+' '.join(map(str, row_rna_ls))

	return savefn

def main():
	sample_RRI_type_dis()
	inter_stat_txt_ls = ['/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-%s-rep-combine/27-DG.inter.stat.txt'%(i) for i in [1,2,3,4,5]]
	for inter_stat_txt in inter_stat_txt_ls:
		RRI_type_dist_conversion(txt=inter_stat_txt)

if __name__ == '__main__':
	main()