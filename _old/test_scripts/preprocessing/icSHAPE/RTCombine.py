"""This script is used to wrap the combine replicates step for all samples.
   It require user to put all replicates in a dictionary and 
   named the two replicates to be combined in the form like [prefix]1 and [prefix]2"""


import subprocess
import itertools
from nested_dict import nested_dict
import sys,os

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

def combine_replicates():
	file_info_dict = file_info("/Share/home/zhangqf7/gongjing/zebrafish/data/RT_raw_new")
	rep1_ls = sorted([i for i in file_info_dict['paths'] if i.endswith('1')])
	rep2_ls = sorted([i for i in file_info_dict['paths'] if i.endswith('2')])
	savefn_dir = '/Share/home/zhangqf7/gongjing/zebrafish/data/RT_combine'
	combine_pl = '/Share/home/zhangqf7/gongjing/zebrafish/script/icSHAPE-master/scripts/combineRTreplicates.pl'
	for rep1, rep2 in zip(rep1_ls, rep2_ls):
		print "combine: %s, %s"%(rep1, rep2)
		sample = rep1.split('/')[-1].replace('_rep1', '')
		savefn = savefn_dir + '/' + sample
		subprocess.call(["perl %s -i %s:%s -o %s"%(combine_pl, rep1, rep2, savefn)], shell=True)


def main():
	combine_replicates()

if __name__ == '__main__':
	main()
