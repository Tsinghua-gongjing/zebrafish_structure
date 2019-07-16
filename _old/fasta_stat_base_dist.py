from pyfasta import Fasta
from nested_dict import nested_dict
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")

def read_fa(fa='/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.fa'):
	fa_dict1 = Fasta(fa, key_fn=lambda key:key.split("\t")[0])
	fa_dict = {i.split()[0].split('.')[0]:j[0:] for i,j in fa_dict1.items()}
	print fa_dict.keys()[0:3]
	return fa_dict

def loadTransGtfBed2(ref_bed='/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.trans.bed2'):
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

def calc_atcg_by_element():
	fa_dict = read_fa()
	trans_dict = loadTransGtfBed2()

	element_atcg_dict = nested_dict(2, int)

	df = pd.read_csv('/Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelationPairwiseNew/2018_03_12_gini/RPKM_combine.merge.txt', header=0, index_col=0, sep='\t')
	print df.index[0:10]
	tx_dict = {i:0 for i in df.index}

	for i,j in trans_dict.items():
		# if not tx_dict.has_key(i):
			# continue
		if int(j['utr_5_end']) > 0:
			utr_5_seq = fa_dict[i][int(j['utr_5_start'])-1:int(j['utr_5_end'])]
			for base in ['A', 'T', 'C', 'G']:
				element_atcg_dict['UTR5'][base] += utr_5_seq.upper().count(base)
		if int(j['cds_end']) > 0:
			cds_seq = fa_dict[i][int(j['cds_start'])-1:int(j['cds_end'])]
			for base in ['A', 'T', 'C', 'G']:
				element_atcg_dict['CDS'][base] += cds_seq.upper().count(base)
		if int(j['utr_3_end']) > 0:
			utr_3_seq = fa_dict[i][int(j['utr_3_start'])-1:int(j['utr_3_end'])]
			for base in ['A', 'T', 'C', 'G']:
				element_atcg_dict['UTR3'][base] += utr_3_seq.upper().count(base)

	print element_atcg_dict

	element_atcg_df = pd.DataFrame.from_dict(element_atcg_dict, orient='index')
	element_atcg_df['total'] = element_atcg_df.sum(axis=1)
	for base in ['A', 'T', 'C', 'G']:
		element_atcg_df['percent(%s)'%(base)] = element_atcg_df[base] / element_atcg_df['total']
	print element_atcg_df

	fig, ax = plt.subplots(figsize=(6,4))
	base_ls = ['percent(%s)'%(base) for base in ['A', 'T', 'C', 'G']]
	element_atcg_df[base_ls].plot(kind='barh', stacked=True)
	plt.tight_layout()
	savefn = '/Share/home/zhangqf7/gongjing/zebrafish/result/base_dist/element_atcg_dist.png'
	plt.savefig(savefn)
	plt.close()


def main():
	calc_atcg_by_element()

if __name__ == '__main__':
	main()