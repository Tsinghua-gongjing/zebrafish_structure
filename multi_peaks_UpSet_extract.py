import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
from nested_dict import nested_dict
import sys,os
import pandas as pd, numpy as np
import itertools

def read_bed(bed):
	print "load: %s"%(bed)
	d = nested_dict(2, list)
	with open(bed, 'r') as BED:
		for line in BED:
			if not line or line.startswith('#'):
				continue
			arr = line.split('\t')
			tx_start_end = '|'.join(arr[1:4])
			sample_ls = arr[6].split('|')
			for sample in sample_ls:
				d[sample][len(sample_ls)].append(tx_start_end)
	d = d.to_dict()
	for sample in sample_ls:
		print sample, len(d[sample][len(sample_ls)]), len(sample_ls)
	return d

def generate_sample(sample_ls=None):
	if sample_ls is None:
		sample_ls = ['egg_cell1', 'cell1_cell4', 'cell4_cell64', 'cell64_sphere', 'sphere_shield']
	all_sample = []
	for i in range(1, len(sample_ls)+1):
		combinations = itertools.combinations(sample_ls, i)
		for c in combinations:
			# s = 'mergePeaks_'+'_'.join([j+'_window-anno.bed' for j in c])
			s = 'new_mergepeaks_d10_'+'_'.join([j+'-var.sorted.merge.bed' for j in c])
			all_sample.append(s)
	print all_sample

	return all_sample

def main():
	# sample_ls = ['egg_cell1', 'cell1_cell4', 'cell4_cell64', 'cell64_sphere', 'sphere_shield']
	sample_ls = ['egg_cell1_egg_cell1', 'cell1_cell4_cell1_cell4', 'cell4_cell64_cell4_cell64', 'cell64_sphere_cell64_sphere', 'sphere_shield_sphere_shield']
	all_sample = generate_sample(sample_ls)
	all_sample_d = nested_dict(2, list)
	# save_dir = '/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/mergepeaks_d10'
	save_dir = '/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/new_mergepeaks_d10'
	for sample in all_sample:
		bed = '%s/%s'%(save_dir, sample)
		d = read_bed(bed)
		for i,j in d.items():
			for m,n in j.items():
				for tx_start_end in n:
					all_sample_d[i][m].append(tx_start_end)
	print len(all_sample_d['egg_cell1/window-anno.bed']), len(all_sample_d['egg_cell1/window-anno.bed'][1]), len(all_sample_d['egg_cell1/window-anno.bed'][4])

	for i,j in all_sample_d.items():
		savefn = '%s/separate/%s.bed'%(save_dir, i.split('/')[0])
		with open(savefn, 'w') as SAVEFN:
			for m,n in j.items():
				for tx_start_end in n:
					print >>SAVEFN, tx_start_end.replace('|', '\t')

	for way in range(1, len(sample_ls)+1):
		way_ls = []
		for i,j in all_sample_d.items():
				for m,n in j.items():
					if m == way:
						for tx_start_end in n:
							way_ls.append(tx_start_end)
		savefn = '%s/separate/way%s.bed'%(save_dir, way)
		with open(savefn, 'w') as SAVEFN:
			for tx_start_end in set(way_ls):
				print >>SAVEFN, tx_start_end.replace('|', '\t')


if __name__ == '__main__':
	main()