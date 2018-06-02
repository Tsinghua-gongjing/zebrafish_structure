import networkx as nx 
from nested_dict import nested_dict
import pandas as pd, numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
import gj
sns.set(style="ticks")
sns.set_context("poster")

def read_rpkm2(rpkm=None):
	if rpkm is None:
		rpkm = '/Share/home/zhangqf7/gongjing/zebrafish/result/RPKMCorrelationPairwiseNew/RPKM_combine.merge.txt'
	df = pd.read_csv(rpkm, header=0, sep='\t', index_col=0)
	rpkm_dict = df.to_dict(orient='index')

	return rpkm_dict

def sample_relation():
	sample_dict = {'shi-zp-1-rep-combine': 'DMSO_egg',
				   'shi-zp-2-rep-combine': 'DMSO_1cell',
				   'shi-zp-3-rep-combine': 'DMSO_4cell',
				   'shi-zp-4-rep-combine': 'DMSO_64cell',
				   'shi-zp-5-rep-combine': 'DMSO_1K'}
	return sample_dict

def loadTransGtfBed2(ref_bed='/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/refseq_ensembl91_merge.tarns.bed.2'):
    H = open(ref_bed)
    line = H.readline()
    trans_dict = nested_dict()
    header_ls = ['tx', 'gene', 'type', 'length', 'utr_5_start', 'utr_5_end', 'cds_start', 'cds_end', 'utr_3_start', 'utr_3_end']
    while line:
        if line.startswith('#'): line = H.readline(); continue
        arr = line.strip().split('\t')
        gene = arr[1].split('=')[0].split()[0]
        for i,j in zip(header_ls, arr):
            trans_dict[arr[0]][i] = j
        line = H.readline()
    H.close()
    print "read: %s, n=%s"%(ref_bed, len(trans_dict))
    return trans_dict.to_dict()

def rpkm_in_RRI_or_not():
	rpkm_dict = read_rpkm2()
	sample_dict = sample_relation()
	trans_dict = loadTransGtfBed2()

	sample_ls = ['shi-zp-1-rep-combine',  'shi-zp-4-rep-combine', ]
	savedir = '/Share/home/zhangqf7/gongjing/zebrafish/result/paris_RRI_abundance'

	df_ls = []
	for sample in sample_ls:
		inter_RRI = '/Share/home/zhangqf7/gongjing/zebrafish/data/paris/%s/27-DG.inter.element.txt'%(sample)
		nx_inter_RRI, df_inter_RRI = read_inter_RRI(inter_RRI, sp=3) # sp: support reads

		RRI_mRNA_ls = [i for i in nx_inter_RRI.nodes() if trans_dict[i]['type'] == 'mRNA']
		RPKM_mRNA_ls = [i for i in rpkm_dict if i.startswith('NM_')]
		non_RRI_mRNA_ls = [i for i in RPKM_mRNA_ls if i not in RRI_mRNA_ls]

		savefn = savedir + '/' + 'RRI_nonRRI_mRNA.%s.%s.txt'%(sample, sample_dict[sample])
		with open(savefn, 'w') as SAVEFN:
			for i in RRI_mRNA_ls:
				if rpkm_dict.has_key(i):
					print >>SAVEFN, '\t'.join(map(str, [i, rpkm_dict[i][sample_dict[sample]], 'RRI mRNA' , sample_dict[sample].split('_')[1]]))
			for i in non_RRI_mRNA_ls:
				if rpkm_dict.has_key(i):
					print >>SAVEFN, '\t'.join(map(str, [i, rpkm_dict[i][sample_dict[sample]], 'non-RRI mRNA', sample_dict[sample].split('_')[1]]))

		df = pd.read_csv(savefn, header=None, sep='\t')
		df.columns = ['tx_id', 'log2(RPKM)', 'mRNA type', 'Stage']
		df_ls.append(df)
	df_all = pd.concat(df_ls)
	print df_all.head(), df_all.shape

	fig,ax = plt.subplots(figsize=(5,4))
	sns.boxplot(x='Stage', y='log2(RPKM)', hue='mRNA type', data=df_all, width=0.4)
	savefn = savedir + '/' + 'RRI_nonRRI_mRNA.abundance.sp3.png'
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.tight_layout()
	plt.savefig(savefn, bbox_inches="tight")
	plt.close()

def read_inter_RRI(inter_RRI=None, sp=3):
	if inter_RRI is None:
		inter_RRI = '/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-2/27-DG.inter.element.txt'
	df_inter_RRI = pd.read_csv(inter_RRI, header=None, sep='\t')
	df_inter_RRI = df_inter_RRI[(df_inter_RRI[13] != 'rRNA') & (df_inter_RRI[14] != 'rRNA') & (df_inter_RRI[9] >= sp)]
	print df_inter_RRI.head()
	nx_inter_RRI = nx.from_pandas_dataframe(df_inter_RRI, 1, 5)
	print "\nread: %s"%(inter_RRI)
	print nx.info(nx_inter_RRI)
	print

	return nx_inter_RRI, df_inter_RRI

def main():
	rpkm_in_RRI_or_not()

if __name__ == '__main__':
	main()