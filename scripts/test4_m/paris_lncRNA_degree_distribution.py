import gj
from nested_dict import nested_dict
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")

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

def read_dg_txt(dg_txt=None, support=3, filter_rRNA=True, only_mRNA_lncRNA=True):
	if dg_txt is None:
		dg_txt = '/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-4/7-DG.txt'
	dg_dict = nested_dict()
	with open(dg_txt, 'r') as DG:
		for line in DG:
			line = line.strip()
			if line.startswith('#'):
				header_ls = line.replace('#', '').split('\t')
				continue
			if not line:
				continue
			arr = line.split('\t')
			if int(arr[9]) < support:
				continue
			if (filter_rRNA and arr[13] == 'rRNA') or (filter_rRNA and arr[14] == 'rRNA'):
				continue
			if only_mRNA_lncRNA and arr[13] not in ['mRNA', 'lncRNA']:
				continue
			if only_mRNA_lncRNA and arr[14] not in ['mRNA', 'lncRNA']:
				continue
			for i,j in zip(header_ls, arr):
				dg_dict[arr[0]][i] = j
			if arr[1] == arr[5]:
				dg_dict[arr[0]]['RRI_type'] = 'intra'
			else:
				dg_dict[arr[0]]['RRI_type'] = 'inter'
	print "DG num: %s, file: %s"%(len(dg_dict), dg_txt)
	return dg_dict.to_dict()

def degree_hist(dg_txt=None):
	if dg_txt is None:
		dg_txt = '/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-5-rep-combine/27-DG.txt'
	dg_dict = read_dg_txt(dg_txt)

	trans_dict = loadTransGtfBed2()

	RRI_dict = nested_dict(2, list)
	for i,j in dg_dict.items():
		if j['lchr'] != j['rchr']:
			RRI_dict[j['RRI_type']][j['lchr']].append(j['rchr'])
			RRI_dict[j['RRI_type']][j['rchr']].append(j['lchr'])
		else:
			RRI_dict[j['RRI_type']][j['lchr']].append(j['rchr'])

	for i in ['inter', 'intra']:
		savefn = dg_txt.replace('.txt', '.%s.degree.txt'%(i))
		degree_ls_ls = [[], [], []]

		with open(savefn, 'w') as SAVEFN:
			for k,v in RRI_dict[i].items():
				print >>SAVEFN, '\t'.join(map(str, [ k, trans_dict[k]['type'], len(v), len(set(v)), ','.join(list(set(v))) ]))
				degree_ls_ls[0].append(len(set(v)))

				if trans_dict[k]['type'] == 'mRNA':
					degree_ls_ls[1].append(len(set(v)))
				if trans_dict[k]['type'] == 'lncRNA':
					degree_ls_ls[2].append(len(set(v)))

		degree_mean_ls = [np.mean(i) for i in degree_ls_ls]
		gj.cumulate_dist_plot(ls_ls=degree_ls_ls,ls_ls_label=['%s, mean=%.2f'%(i,j) for i,j in zip(['all', 'mRNA', 'lncRNA'], degree_mean_ls)], bins=40,title='degree distribution',ax=None,savefn=savefn.replace('.txt', '.pdf'),xlabel='log2(# of interacting partners)',ylabel=None,add_vline=None,add_hline=None,log2transform=1,xlim=None,ylim=None)

def main():
	dg_txt_ls = ['/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-%s-rep-combine/27-DG.txt'%(i) for i in [1,3,4,5]]
	for dg_txt in dg_txt_ls:
		degree_hist(dg_txt)

if __name__ == '__main__':
	main()
