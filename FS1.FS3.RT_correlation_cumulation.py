import sys
# sys.path.append('/Share/home/zhangqf7/gongjing/zebrafish/script')
import matplotlib as mpl
mpl.use('Agg')
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import gj
from glob import glob



ls_ls, ls_ls_label = [], []
file_ls = sys.argv[1:-1]
savefn = sys.argv[-1]
print "input: ", file_ls
print "output: ", savefn
for sample_fn in file_ls:
    df = pd.read_csv(sample_fn, header=None, sep='\t',  keep_default_na=False, na_values=['n/a'])
    df.dropna(axis=0, how='any', inplace=True)
    ls_ls.append(list(df[4]))
    ls_ls_label.append(sample_fn.split('/')[-1])
#    savefn = './test.png'
gj.cumulate_dist_plot(ls_ls=ls_ls,ls_ls_label=ls_ls_label,bins=40,title=None,ax=None,savefn=savefn,xlabel=None,ylabel='',add_vline=[0.6, 0.7],add_hline=[0,0.3, 0.4,1],log2transform=0, xlim=[-0.05,1.05], ylim=[-0.05,1.05])
