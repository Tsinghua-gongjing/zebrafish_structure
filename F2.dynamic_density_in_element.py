import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
from nested_dict import nested_dict
import sys,os
import pandas as pd, numpy as np
from pyfasta import Fasta

def readIc(ic_file):
	"Read icSHAPE Reactivity File"
	icDict = dict()
	IN = open(ic_file)
	line = IN.readline()
	while line:
		arr = line.strip().split()
		icDict[arr[0]] = arr[3:]
		line = IN.readline()
	IN.close()
	return icDict

def read_fa(fa='/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.fa'):
    fa_dict1 = Fasta(fa, key_fn=lambda key: key.split("\t")[0])
    fa_dict = {i.split()[0]: j[0:] for i, j in fa_dict1.items()}
    print fa_dict.keys()[0:3]
    return fa_dict

def loadTransGtfBed2(ref_bed='/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/ref_GRCz10_top_level.trans.bed.2'):
	H = open(ref_bed)
	line = H.readline()
	trans_dict = nested_dict()
	header_ls = ['tx', 'gene', 'type', 'length', 'utr_5_start', 'utr_5_end', 'cds_start', 'cds_end', 'utr_3_start', 'utr_3_end']
	while line:
		if line.startswith('#'): line = H.readline(); continue
		arr = line.strip().split()
		for i,j in zip(header_ls, arr):
			trans_dict[arr[0]][i] = j
		line = H.readline()
	H.close()
	print "read: %s, n=%s"%(ref_bed, len(trans_dict))
	return trans_dict.to_dict()

def list_split_equal(list_of_values=None, split_num=10, mode='gini', null_pct=1):
	array_parts = np.array_split(list_of_values, split_num)
	array_parts = [[i[0], i[-1]] for i in array_parts]
	return array_parts

def ic_density(bed=None, cut_num_ls=None, savefn=None, split_overlap_ratio_min=0.5, sample='cell1_cell4'):
	if bed is None:
		bed = '/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/%s/window-anno.bed'%(sample)
	if cut_num_ls is None:
		cut_num_ls = [20, 200, 80]
	if savefn is None:
		savefn = '/Share/home/zhangqf7/gongjing/zebrafish/result/icshape_signal_mean/sample_%s_dynamic_density.txt'%(sample)

	trans_dict = loadTransGtfBed2('/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.trans.bed2')

	bed_region_dict = nested_dict(1, list)
	with open(bed, 'r') as BED:
		for line in BED:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			arr = line.split('\t')
			if [int(arr[1]), int(arr[2])] not in bed_region_dict[arr[0]]:
			# if [arr[1], arr[2]] not in bed_region_dict[arr[0]]:  # [arr[1], arr[2]] is str, while [int(arr[1]), int(arr[2])] is int
				bed_region_dict[arr[0]].append([int(arr[1]), int(arr[2])])
	print bed_region_dict['NM_205538']

	anno = bed.replace('.bed', '.element.txt')
	ANNO = open(anno, 'w')

	fa_dict = read_fa()
	seq_savefn = savefn.replace('.txt', '.seq.txt')
	SEQ = open(seq_savefn, 'w')

	with open(savefn, 'w') as SAVEFN:
		for tx,j in bed_region_dict.items():
			if not trans_dict.has_key(tx):
				continue
			tx_element_count = [0]* sum(cut_num_ls)
			utr_5_start, utr_5_end, cds_start, cds_end, utr_3_start, utr_3_end = [int(trans_dict[tx][g]) for g in ['utr_5_start', 'utr_5_end', 'cds_start', 'cds_end', 'utr_3_start', 'utr_3_end']]
			if utr_5_end < cut_num_ls[0]:
				continue
			if cds_end - cds_start + 1 < cut_num_ls[1]:
				continue
			if utr_3_end - utr_3_start + 1 < cut_num_ls[2]:
				continue
			utr_5_split = list_split_equal(xrange(utr_5_start, utr_5_end+1), cut_num_ls[0])
			cds_split = list_split_equal(xrange(cds_start, cds_end+1), cut_num_ls[1])
			utr_3_split = list_split_equal(xrange(utr_3_start, utr_3_end+1), cut_num_ls[2])

			all_split = utr_5_split + cds_split + utr_3_split # 1-based

			for (j_start, j_end) in j:
				for n,(split_start, split_end) in enumerate(all_split):
					if max(j_start, split_start) < min(j_end, split_end):
						overlap_len = min(j_end, split_end) - max(j_start, split_start) + 1
						split_overlap_ratio = overlap_len / float(split_end-split_start+1)
						if split_overlap_ratio >= split_overlap_ratio_min:
							tx_element_count[n] += 1
							print >>ANNO, '\t'.join(map(str, [tx, j_start, j_end, split_start, split_end, n, overlap_len, split_overlap_ratio]))

			print >>SAVEFN, '\t'.join(map(str, [tx]+tx_element_count))

			tx_a_content_cout = [0]* sum(cut_num_ls)
			tx_t_content_cout = [0]* sum(cut_num_ls)
			tx_c_content_cout = [0]* sum(cut_num_ls)
			tx_g_content_cout = [0]* sum(cut_num_ls)
			for n,(split_start, split_end) in enumerate(all_split):
				a = fa_dict[tx][split_start-1:split_end].upper().count('A') / float(len(fa_dict[tx][split_start-1:split_end].upper()))
				t = fa_dict[tx][split_start-1:split_end].upper().count('T') / float(len(fa_dict[tx][split_start-1:split_end].upper()))
				c = fa_dict[tx][split_start-1:split_end].upper().count('C') / float(len(fa_dict[tx][split_start-1:split_end].upper()))
				g = fa_dict[tx][split_start-1:split_end].upper().count('G') / float(len(fa_dict[tx][split_start-1:split_end].upper()))
				tx_a_content_cout[n] = a
				tx_t_content_cout[n] = t
				tx_c_content_cout[n] = c
				tx_g_content_cout[n] = g
			print >>SEQ, '\t'.join(map(str, [tx, 'A']+tx_a_content_cout)) 
			print >>SEQ, '\t'.join(map(str, [tx, 'T']+tx_t_content_cout)) 
			print >>SEQ, '\t'.join(map(str, [tx, 'C']+tx_c_content_cout)) 
			print >>SEQ, '\t'.join(map(str, [tx, 'G']+tx_g_content_cout)) 

	ANNO.close()
	SEQ.close()

	return savefn

def ic_density_plot(txt=None):
	if txt is None:
		txt = '/Share/home/zhangqf7/gongjing/zebrafish/result/icshape_signal_mean/all_sample_dynamic_density.txt'
	df = pd.read_csv(txt, header=None, sep='\t', index_col=0)
	df_stat = df.sum()
	print df_stat
	df_stat = df_stat / float(df_stat.sum())
	fig,ax=plt.subplots()
	df_stat.plot(kind='line', ax=ax)
	cut_num_ls = [20, 200, 80]
	plt.axvline(x=cut_num_ls[0]+0.5,linestyle='--',color='black')
	plt.axvline(x=cut_num_ls[0]+cut_num_ls[1]+0.5,linestyle='--',color='black')
	plt.savefig(txt.replace('.txt', '.pdf'))
	plt.tight_layout()
	plt.close()

def ic_density_plot_all(sample_ls=None, savefn=None):
	if sample_ls is None:
		sample_ls = ['egg_cell1', 'cell1_cell4', 'cell4_cell64', 'cell64_sphere', 'sphere_shield']
		# sample_ls = ['way1', 'way2', 'way3', 'way4', 'way5']
	
	fig,ax=plt.subplots()
	for sample in sample_ls:
		# df = pd.read_csv('/Share/home/zhangqf7/gongjing/zebrafish/result/icshape_signal_mean/sample_%s_dynamic_density.txt'%(sample), header=None, sep='\t', index_col=0)
		df = pd.read_csv('/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/mergepeaks_d10/separate/%s.bed.txt'%(sample), header=None, sep='\t', index_col=0)
		df_stat = df.sum()
		print df_stat
		df_stat = df_stat / float(df_stat.sum())
		df_stat.plot(kind='line', ax=ax, label='%s'%(sample))
	plt.axvline(x=20.5,linestyle='--',color='black')
	plt.axvline(x=220.5,linestyle='--',color='black')
	plt.legend()
	plt.tight_layout()
	plt.savefig(savefn)
	plt.close()

def ic_seq_plot(txt=None):
	if txt is None:
		txt = '/Share/home/zhangqf7/gongjing/zebrafish/result/icshape_signal_mean/all_sample_dynamic_density.seq.txt'
	cut_num_ls = [10, 100, 40]
	df = pd.read_csv(txt, header=None, sep='\t', index_col=0)
	df.columns = ['base'] + range(1, sum(cut_num_ls)+1)
	print df.head()

	fig,ax=plt.subplots()
	for b in ['A', 'T', 'C', 'G']:
		df_base = df[df['base']==b]
		print df_base.head()
		df_base.drop(['base'], inplace=True, axis=1)
		df_stat = df_base.mean()
		df_stat.plot(kind='line', ax=ax, label=b)

	plt.axvline(x=cut_num_ls[0]+0.5,linestyle='--',color='black')
	plt.axvline(x=cut_num_ls[0]+cut_num_ls[1]+0.5,linestyle='--',color='black')
	plt.legend()
	plt.savefig(txt.replace('.txt', '.png'))
	plt.tight_layout()
	plt.close()

def main():
	# savefn = ic_density(sample='egg_cell1')
	savefn = ic_density(bed='/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/dynamic_region/way2345/way2345.bed', 
			   savefn='/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/dynamic_region/way2345/way2345.bed.txt')
	ic_density_plot(savefn)

	savefn = ic_density(bed='/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/dynamic_region/way1/way1.bed', 
			   savefn='/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/dynamic_region/way1/way1.bed.txt')
	ic_density_plot(savefn)
	# ic_seq_plot('/Share/home/zhangqf7/gongjing/zebrafish/result/icshape_signal_mean/sample_egg_cell1_dynamic_density.seq.txt')

	# ic_density_plot_all()
	# ic_density_plot_all(sample_ls=['way1', 'way2', 'way3', 'way4', 'way5'], savefn='/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/mergepeaks_d10/separate/way_all.pdf')
	# ic_density_plot_all(savefn='/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/mergepeaks_d10/separate/sample_individual_all.pdf')

if __name__ == '__main__':
	main()