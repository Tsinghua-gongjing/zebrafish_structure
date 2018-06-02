import networkx as nx 
import pandas as pd
from nested_dict import nested_dict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
sns.set(style="ticks")
sns.set_context("poster")
import venn_all
import os
import numpy as np

def read_inter_RRI(inter_RRI=None, filter_rRNA=False, support_read_num=2, only_mRNA_lncRNA=False):
	if inter_RRI is None:
		inter_RRI = '/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-2/27-DG.inter.element.txt'
	df_inter_RRI = pd.read_csv(inter_RRI, header=None, sep='\t')
	if filter_rRNA:
		df_inter_RRI = df_inter_RRI[(df_inter_RRI[13] != 'rRNA') & (df_inter_RRI[14] != 'rRNA')]
	if only_mRNA_lncRNA:
		only_mRNA_lncRNA_index = (df_inter_RRI[13].isin(['mRNA', 'lncRNA'])) & (df_inter_RRI[14].isin(['mRNA', 'lncRNA']))
		df_inter_RRI = df_inter_RRI[only_mRNA_lncRNA_index]
	header_ls = ['Group', 'lchr', 'lstrand', 'lstart', 'lend', 'rchr', 'rstrand',
				'rstart', 'rend', 'support', 'lcount', 'rcount', 'score', 'ltype', 'rtype', 'RRI_type', 'lcontext', 'rcontext']
	df_inter_RRI.columns = header_ls
	df_inter_RRI = df_inter_RRI[df_inter_RRI['support'] >= support_read_num]
	print df_inter_RRI.head()
	nx_inter_RRI = nx.from_pandas_dataframe(df_inter_RRI, 'lchr', 'rchr', edge_attr=['support', 'Group'])
	print "\nread: %s"%(inter_RRI)
	print nx.info(nx_inter_RRI)
	print

	return nx_inter_RRI, df_inter_RRI

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

def read_maternal(maternal_list_file=None):
	if maternal_list_file is None:
		matenal = '/Share/home/zhangqf7/gongjing/zebrafish/data/maternal_gene/maternal-decay.txt'
	else:
		matenal = maternal_list_file
	matenal_gene_dict = nested_dict(1, int)
	with open(matenal, 'r') as IN:
		for line in IN:
			line = line.strip()
			if not line or line.startswith('value'):
				continue
			arr = line.split('\t')
			matenal_gene_dict[arr[0]] += 1
	return matenal_gene_dict.to_dict()


def merge_RRI_new(stage_ls, sample_ls):
	inter_RRI_ls = ['/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-%s-rep-combine/27-DG.inter.element.txt'%(i) for i in stage_ls]
	nx_inter_RRI_network_ls = []
	for inter_RRI,sample in zip(inter_RRI_ls, sample_ls):
		nx_inter_RRI, df_inter_RRI = read_inter_RRI(inter_RRI, only_mRNA_lncRNA=True, support_read_num=3)
		nx_inter_RRI_network_ls.append(nx_inter_RRI)
	nx_inter_RRI_merge = nx.compose_all(nx_inter_RRI_network_ls)

	print nx.info(nx_inter_RRI_merge)

	trans_dict = loadTransGtfBed2()
	maternal_decay_file = '/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/decay.txt'
	maternal_decay_dict = read_maternal(maternal_list_file=maternal_decay_file)
	maternal_stable_file = '/Share/home/zhangqf7/gongjing/zebrafish/result/zhangting_provided/supplementary/stable.txt'
	maternal_stable_dict = read_maternal(maternal_list_file=maternal_stable_file)

	def _decayorstable(tx_id):
		
		if tx_id in maternal_decay_dict:
			if tx_id in maternal_stable_dict:
				return 'mystic'
			else:
				return 'decay'
		elif tx_id in maternal_stable_dict:
			return 'stable'
		else:
			return 'undetermined'
	out_dir = '/Share/home/zhangqf7/gongjing/zebrafish/result/network/merge/_select_example_for_fig4'
	savefn = os.path.join(out_dir, 'nx_merge_%s.mrna_lncrna.nodes.new.decay.txt' % "".join(map(str,stage_ls)))
	savefn2 = os.path.join(out_dir, 'nx_merge_%s.mrna_lncrna.edges.new.decay.txt' % "".join(map(str,stage_ls)))
	SAVEFN2 = open(savefn2, 'w')
	SAVEFN = open(savefn, 'w')
	header_ls = ['tx_id', 'gene_name', 'type', 'Maternal_Decay'] + sample_ls
	print >>SAVEFN, '\t'.join(header_ls)
	for node in nx_inter_RRI_merge.nodes():
		node_status = ['YES' if n.has_node(node) else 'NO' for n in nx_inter_RRI_network_ls]
		maternalorNot = _decayorstable(node)
		print >>SAVEFN, '\t'.join([node, trans_dict[node]['gene'].split('=')[0], trans_dict[node]['type'], maternalorNot] + node_status)
	SAVEFN.close()

	header_ls2 = ['source', 'target'] + sample_ls + [ "%s(read)" % i for i in sample_ls] + ["%s(group)" % i for i in sample_ls] + ['type_source', 'type_target', 'name_source', 'name_target', 
				'source_Maternal_Decay', 'target_Maternal_Decay']
	print >>SAVEFN2, '\t'.join(header_ls2)
	for i in nx_inter_RRI_merge.edges():
		edge_status = ['YES' if n.has_edge(i[0], i[1]) else 'NO' for n in nx_inter_RRI_network_ls]
		edge_support_read = [n[i[0]][i[1]]['support'] if n.has_edge(i[0], i[1]) else 0 for n in nx_inter_RRI_network_ls]
		edge_group = [n[i[0]][i[1]]['Group'] if n.has_edge(i[0], i[1]) else '.' for n in nx_inter_RRI_network_ls]
		print >> SAVEFN2, '\t'.join([i[0], i[1]] + edge_status + map(str, edge_support_read) + map(str, edge_group) + [trans_dict[i[0]]['type'], trans_dict[i[1]]['type']]
									+ [trans_dict[i[0]]['gene'].split('=')[0], trans_dict[i[1]]['gene'].split('=')[0]] + [_decayorstable(i[0]), _decayorstable(i[1])])
	SAVEFN2.close()


def main():
    merge_RRI_new([1,4], ['egg','64cell'])


if __name__ == '__main__':
	main()
