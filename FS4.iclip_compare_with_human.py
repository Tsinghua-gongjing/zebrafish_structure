import pandas as pd 
from nested_dict import nested_dict
import subprocess

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
plt.rcParams["font.family"] = "helvetica"

import venn_all

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

def read_ortholog(fn='/Share2/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/transcript_human_ortholog.txt'):
	ortholog_dict = nested_dict()
	df = pd.read_csv(fn, header=0, sep='\t')
	for i,j in zip(df['Gene name'], df['Human gene name']):
		ortholog_dict[i] = j
	return ortholog_dict.to_dict()

def read_ortholog2(fn='/Share2/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/transcript_human_ortholog.txt'):
	ortholog_dict = nested_dict()
	df = pd.read_csv(fn, header=0, sep='\t')
	for i,j in zip(df['Gene name'], df['Human gene name']):
		ortholog_dict[j] = i
	return ortholog_dict.to_dict()

def read_bed(bed, trans_dict, ortholog_dict):
	df = pd.read_csv(bed, header=None, sep='\t')
	NM_ls = list(df[0])
	gene_ls = [trans_dict[i]['gene'].split('=')[0] if trans_dict.has_key(i) else i for i in NM_ls]
	orth_gene_ls = [ortholog_dict[i] if ortholog_dict.has_key(i) else i for i in gene_ls]
	return orth_gene_ls

def read_human_hur():
	fn1 = '/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/human_HUR/1-s2.0-S1097276511004229-mmc3.xls'
	df1 = pd.read_excel(fn1, sheet_name='consensus')
	gene_ls1 = df1['Gene']
	print "table1: %s"%(len(gene_ls1))

	fn2 = '/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/human_HUR/molcel_3918_mmc3.csv'
	df2 = pd.read_csv(fn2, header=0, sep=',')
	gene_ls2 = df2[df2['Total_Binding_Sites']>=1]['Symbol']
	print "table2: %s"%(len(gene_ls2))

	return gene_ls1, gene_ls2

trans_dict = loadTransGtfBed2()
ortholog_dict = read_ortholog()

human1, human2 = read_human_hur()
h4 = read_bed('/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_all_rep12.c6.utr3.bed', trans_dict, ortholog_dict)
h6 = read_bed('/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h6_all_rep12.c6.utr3.bed', trans_dict, ortholog_dict)

labels = venn_all.get_labels(data=[set(human1), set(human2), set(h4), set(h6)])
labels_ls = ['human1', 'human2', 'H4', 'H6']
venn_all.venn4(labels=labels, names=labels_ls, figsize=(6,6))
plt.tight_layout()
savefn = '/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_h6_rep12.c6.humanHuR.pdf'
plt.savefig(savefn)
plt.close()


ortholog_dict2 = read_ortholog2()
common = set(human1) & set(human2) & set(h4) & set(h6)
with open(savefn.replace('.pdf', '.txt'), 'w') as FN:
	for i in common:
		print >>FN, '\t'.join([i, ortholog_dict2[i]])

labels = venn_all.get_labels(data=[set(human1), set(h4)])
labels_ls = ['human1', 'H4']
venn_all.venn2(labels=labels, names=labels_ls, figsize=(6,6))
plt.tight_layout()
savefn = '/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_h6_rep12.c6.human1HuR.pdf'
plt.savefig(savefn)
plt.close()

labels = venn_all.get_labels(data=[set(human2), set(h4)])
labels_ls = ['human2', 'H4']
venn_all.venn2(labels=labels, names=labels_ls, figsize=(6,6))
plt.tight_layout()
savefn = '/Share2/home/zhangqf7/gongjing/zebrafish/result/iCLIP/h4_h6_rep12.c6.human2HuR.pdf'
plt.savefig(savefn)
plt.close()

		