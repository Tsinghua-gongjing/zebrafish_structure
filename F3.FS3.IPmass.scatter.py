from pandas.io.parsers import read_csv
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
import statsmodels.stats.multitest as multi
from adjustText import adjust_text

def parse_IPMASS(t=None, mode='log2'):
    if t is None:
        t = 'T1'
    f = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/IP-mass/%s_P_N.xlsx'%(t)

    df = pd.read_excel(f)
    print "read: %s, num=%s"%(f, df.shape[0])

#     cols_keep = ['Gene names', #'Q-value', 'Score', 
#                  'LFQ intensity %s_N1'%(t), 'LFQ intensity %s_N2'%(t), 'LFQ intensity %s_N3'%(t), 'LFQ intensity %s_N4'%(t),
#                  'LFQ intensity %s_P1'%(t), 'LFQ intensity %s_P2'%(t), 'LFQ intensity %s_P3'%(t), 'LFQ intensity %s_P4'%(t),]
#     rep_N = ['LFQ intensity %s_N1'%(t), 'LFQ intensity %s_N2'%(t), 'LFQ intensity %s_N3'%(t), 'LFQ intensity %s_N4'%(t),]
#     rep_P = ['LFQ intensity %s_P1'%(t), 'LFQ intensity %s_P2'%(t), 'LFQ intensity %s_P3'%(t), 'LFQ intensity %s_P4'%(t),]
    rep_N = ['LFQ intensity %s_N1'%(t), 'LFQ intensity %s_N2'%(t)]
    rep_P = ['LFQ intensity %s_P1'%(t), 'LFQ intensity %s_P2'%(t)]
    cols_keep = ['Gene names'] + rep_N + rep_P

    # print df.head()

    # log2 first
    if mode == 'log2':
        for i in rep_N + rep_P:
            log2_ls = []
            for v in list(df[i]):
                if float(v) == 0:
                    log2_ls.append(0.001)
                else:
                    log2_ls.append(np.log2(float(v)))
            df['log2(%s)'%(i)] = log2_ls

        df['sum(N)'] = df.loc[:, ['log2(%s)'%(i) for i in rep_N]].sum(axis=1)
        df['sum(P)'] = df.loc[:, ['log2(%s)'%(i) for i in rep_P]].sum(axis=1)
        df['mean(N)'] = df['sum(N)'] / 4.0
        df['mean(P)'] = df['sum(P)'] / 4.0
        df['sum(P)-sum(N)'] = df['sum(P)'] - df['sum(N)']
        df['mean(P)-mean(N)'] = df['mean(P)'] - df['mean(N)']

    #     print df.head()

        rep_N_log2_ls = ['log2(%s)'%(i) for i in rep_N]
        rep_P_log2_ls = ['log2(%s)'%(i) for i in rep_P]

        pvalue = []
        for index, row in df.iterrows():
            rep_N_val_ls = [row[i] for i in rep_N_log2_ls]
            rep_P_val_ls = [row[i] for i in rep_P_log2_ls]
            s,p = stats.ttest_ind(rep_P_val_ls, rep_N_val_ls)
            if np.isnan(p):
                p = 1
            pvalue.append(p)
    else:
        df['sum(N)'] = df.loc[:, ['%s'%(i) for i in rep_N]].sum(axis=1)
        df['sum(P)'] = df.loc[:, ['%s'%(i) for i in rep_P]].sum(axis=1)
        df['mean(N)'] = df['sum(N)'] / 4.0
        df['mean(P)'] = df['sum(P)'] / 4.0
        df['mean(N)'] = [1 if i==0 else i for i in df['mean(N)']]
        df['mean(P)'] = [1 if i==0 else i for i in df['mean(P)']]
        df['log2(mean(N))'] = np.log2(df['mean(N)'])
        df['log2(mean(P))'] = np.log2(df['mean(P)'])
        df['sum(P)-sum(N)'] = df['sum(P)'] - df['sum(N)']
        df['mean(P)-mean(N)'] = df['mean(P)'] - df['mean(N)']
        df['log2(mean(P)/mean(N))'] = df['log2(mean(P))'] - df['log2(mean(N))']
        
        rep_N_log2_ls = ['log2(%s)'%(i) for i in rep_N]
        rep_P_log2_ls = ['log2(%s)'%(i) for i in rep_P]
        
        pvalue = []
        for index, row in df.iterrows():
            rep_N_val_ls = [row[i] for i in rep_N]
            rep_P_val_ls = [row[i] for i in rep_P]
            s,p = stats.ttest_ind(rep_P_val_ls, rep_N_val_ls)
            if np.isnan(p):
                p = 1
            pvalue.append(p)
            
        

    # print pvalue
    qvalue = multi.multipletests(pvalue)
    # print qvalue

    df['pvalue'] = pvalue
    df['qvalue'] = qvalue[1]
    df['-log10(qvalue)'] = -np.log10(df['qvalue'])
    df['-log10(pvalue)'] = -np.log10(df['pvalue'])
    cols_calc = ['sum(N)', 'sum(P)', 'mean(N)', 'mean(P)', 'log2(mean(N))', 'log2(mean(P))', 
                 'sum(P)-sum(N)', 'mean(P)-mean(N)', 'log2(mean(P)/mean(N))',
                'pvalue', 'qvalue', '-log10(pvalue)', '-log10(qvalue)']
    df = df[cols_keep+cols_calc]
    df.to_excel('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/IP-mass/%s_enrich_table.xlsx'%(t), header=True, index=False)
    df.head()
    
    fig,ax=plt.subplots(figsize=(6, 6))
    x_col = 'log2(mean(P)/mean(N))'
    x_col = 'log2(mean(N))'
    y_col = '-log10(pvalue)'
    y_col = 'log2(mean(P))'
    df.plot(kind='scatter', x=x_col, y=y_col, ax=ax)
    ratio_max = max(df[x_col])
#     plt.axvline(x=0, ymin=0, ymax=1, ls='--', color='grey')
#     plt.axhline(y=-np.log10(0.05), xmin=0, xmax=1, ls='--', color='grey')
#     ax.set_xlim(-ratio_max-1, ratio_max+1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.title('time: %s (n=%s)'%(t, df.shape[0]))
    savefn = '/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/IP-mass/%s_enrich_pvalue.pdf'%(t)
    
    df['Gene Names'] = [i.split(';')[0] for i in df['Gene names']]
    texts = []
    for x,y,t in zip(df[x_col], df[y_col], df['Gene Names']):
        if y > -np.log10(0.05) and t == 'elavl1':
#             ax.annotate(t, (x, y), fs=3)
            texts.append(plt.text(x, y, t, fontsize=12))
#     plt.tight_layout()
    adjust_text(texts, only_move={'text': 'x'})
    
    plt.tight_layout()
    plt.savefig(savefn)
    plt.close()
    
    return df, rep_N_log2_ls, rep_P_log2_ls, df[cols_keep+cols_calc]

df1, rep_N_log2_ls1, rep_P_log2_ls1, df1_save = parse_IPMASS(t='T1', mode='normal')
df2, rep_N_log2_ls2, rep_P_log2_ls2, df2_save = parse_IPMASS(t='T2', mode='normal')

def plot_pair(df1, df2, label1, label2, RBP_type_dict):
    matplotlib_venn.venn2([set(list(df1['Gene Names'])), set(list(df2['Gene Names']))], set_labels=[label1, label2])
    df = pd.merge(df1, df2, on='Gene Names', how='outer')
    df.fillna(0, inplace=True)
    df = df[df['Gene Names'].isin(RBP_type_dict.keys()+['pl10'])]
    print df.shape, df1.shape, df2.shape
    print df.columns
#     print df
    
    df = df[(df['log2(mean(P)/mean(N))_x']>=0) & (df['log2(mean(P)/mean(N))_y']>=0)]
    df['RBP_type'] = [RBP_type_dict[i] if RBP_type_dict.has_key(i) else 'other' for i in df['Gene Names']]
    print df.shape
    with open('./%s_%s.genename.txt'%(label1, label2), 'w') as TXT:
        for i in df['Gene Names']:
            print >>TXT, i
    
    fig, ax = plt.subplots(figsize=(14,14))
#     df.plot(kind='scatter', x='log2(mean(P)/mean(N))_x', y='log2(mean(P)/mean(N))_y', ax=ax)
    groups = df.groupby('RBP_type')
    for name, group in groups:
        ax.plot(group['log2(mean(P)/mean(N))_x'], group['log2(mean(P)/mean(N))_y'], marker='o', linestyle='', ms=12, label=name)
    
    ax.legend()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_xlabel('%s, log2(P/N)'%(label1))
    ax.set_ylabel('%s, log2(P/N)'%(label2))
#     ax.set_xlim(0,)
#     ax.set_ylim(0,)

#     df['Gene Names'] = [i.split(';')[0] for i in df['Gene names_x']]
    texts = []
#     for x,y,t in zip(df['log2(mean(P)/mean(N))_x'], df['log2(mean(P)/mean(N))_y'], df['Gene Names']):
#             texts.append(plt.text(x, y, t, fontsize=12))
#     adjust_text(texts, only_move={'text': 'xy'})
    plt.savefig('./%s_%s.enrich_not.pdf'%(label1, label2))

    return df

df = plot_pair(df1, df2, 'T1', 'T2', RBP_type_dict)
df.to_csv('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/IP-mass/T1_T2.txt', header=True, index=True, sep='\t')




# FS.DE
# df2['RBP_type'] = [RBP_type_dict[i] if RBP_type_dict.has_key(i) else 'other' for i in df2['Gene Names']]

# fig, ax = plt.subplots(figsize=(8,14))
# groups = df2.groupby('RBP_type')
# for name, group in groups:
#     ax.plot(group['log2(mean(P)/mean(N))'], group['-log10(qvalue)'], marker='o', linestyle='', ms=12, label=name)
    
# ax.legend()
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# ax.set_xlabel('%s, log2(P/N)'%('P/N'))
# ax.set_ylabel('%s, log2(P/N)'%('-log10(pvalue)'))
# #     ax.set_xlim(0,)
# ax.set_ylim(-0.1,)
# plt.savefig('/Share2/home/zhangqf7/gongjing/zebrafish/script/zebrafish_structure/data/IP-mass/%s.enrich_pvalue.pdf'%('T2'))