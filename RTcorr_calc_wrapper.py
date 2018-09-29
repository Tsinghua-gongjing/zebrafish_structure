from nested_dict import nested_dict
import sys,os
import pandas as pd, numpy as np
import subprocess


def file_info(file_dir=None, result_dir=None):
    if file_dir is None:
        file_dir = '/Share/home/zhangqf7/gongjing/zebrafish/data/RT'
    if result_dir is None:
        result_dir = '/Share/home/zhangqf7/gongjing/zebrafish/result/RTCorrelationAccmulate'
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

def icshape_correlationRT(selection_strategy='useselfbackground', threshold_coverage=2):
    file_info_dict = file_info(file_dir='/Share/home/zhangqf7/gongjing/zebrafish/data/RT_raw_new', result_dir='/Share/home/zhangqf7/gongjing/zebrafish/result/RTCorrelationAccmulate')
    threshold_window_ls = [200]
    icshape_correlatioRT_pl = '/Share/home/zhangqf7/gongjing/zebrafish/script/icSHAPE-master/scripts/correlationRT_test.pl'
    RT1_ls = sorted([i for i in file_info_dict['paths'] if i.endswith('1')])[0:]
    RT2_ls = sorted([i for i in file_info_dict['paths'] if i.endswith('2')])[0:]
    for n,(input_signal_file1,input_signal_file2) in enumerate(zip(RT1_ls, RT2_ls)):
        for threshold_window in threshold_window_ls:
            sample1 = input_signal_file1.split('/')[-1]
            sample2 = input_signal_file2.split('/')[-1]
            output = file_info_dict['result_dir'] + '/' + sample1+'-'+sample2+'-T%st%s'%(threshold_coverage, threshold_window)
            if n%2 == 0 or n%3 == 0:
            	# print "%s -1 %s -2 %s -s %s -f %s -T %s -t %s > %s &"%(icshape_correlatioRT_pl, input_signal_file1, input_signal_file2, selection_strategy, 'rt', threshold_coverage, threshold_window, output)
                subprocess.call(["%s -1 %s -2 %s -s %s -f %s -T %s -t %s > %s &"%(icshape_correlatioRT_pl, input_signal_file1, input_signal_file2, selection_strategy, 'rt', threshold_coverage, threshold_window, output)], shell=True)
            else:
            	# print "%s -1 %s -2 %s -s %s -f %s -T %s -t %s > %s"%(icshape_correlatioRT_pl, input_signal_file1, input_signal_file2, selection_strategy, 'rt', threshold_coverage, threshold_window, output)
                subprocess.call(["%s -1 %s -2 %s -s %s -f %s -T %s -t %s > %s"%(icshape_correlatioRT_pl, input_signal_file1, input_signal_file2, selection_strategy, 'rt', threshold_coverage, threshold_window, output)], shell=True)

if __name__ == '__main__':
	icshape_correlationRT()