import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Helvetica"
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
from nested_dict import nested_dict
import sys,os
import pandas as pd, numpy as np
from scipy import stats

# h4_full_12_all, count>=6
fig, ax = plt.subplots(1,2,figsize=(12,6))

sizes = [532, 4972, 16188]
labels = ['UTR5', 'CDS', 'UTR3']
total = sum(sizes)
ax[0].pie(sizes, labels=labels, 
        autopct=lambda(p): '{:.0f}'.format(p * total / 100), shadow=False, startangle=140)
ax[0].axis('equal')

sizes = [796, 6972, 17700]
labels = ['UTR5', 'CDS', 'UTR3']
total = sum(sizes)
ax[1].pie(sizes, labels=labels, 
        autopct=lambda(p): '{:.0f}'.format(p * total / 100), shadow=False, startangle=140)
ax[1].axis('equal')

plt.savefig('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/FS4.iclip_tag_dist_h4_h6_rep12.c6.pdf')
plt.show()


# # h6_full_12_all, count>=6
# fig, ax = plt.subplots(figsize=(6,6))

# sizes = [796, 6972, 17700]
# labels = ['UTR5', 'CDS', 'UTR3']
# total = sum(sizes)
# ax.pie(sizes, labels=labels, 
#         autopct=lambda(p): '{:.0f}'.format(p * total / 100), shadow=False, startangle=140)
# ax.axis('equal')


# plt.savefig('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/results/iclip_tag_dist_h6_rep12.c6.pdf')
# plt.show()