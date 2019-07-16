def parse_IPMASS(t=None, mode='log2'):
    if t is None:
        t = 'T1'
    f = './%s_P_N.xlsx'%(t)

    df = pd.read_excel(f)
    print "read: %s, num=%s"%(f, df.shape[0])

    cols_keep = ['Gene names', #'Q-value', 'Score', 
                 'LFQ intensity %s_N1'%(t), 'LFQ intensity %s_N2'%(t), 'LFQ intensity %s_N3'%(t), 'LFQ intensity %s_N4'%(t),
                 'LFQ intensity %s_P1'%(t), 'LFQ intensity %s_P2'%(t), 'LFQ intensity %s_P3'%(t), 'LFQ intensity %s_P4'%(t),]
    rep_N = ['LFQ intensity %s_N1'%(t), 'LFQ intensity %s_N2'%(t), 'LFQ intensity %s_N3'%(t), 'LFQ intensity %s_N4'%(t),]
    rep_P = ['LFQ intensity %s_P1'%(t), 'LFQ intensity %s_P2'%(t), 'LFQ intensity %s_P3'%(t), 'LFQ intensity %s_P4'%(t),]

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
    df.to_excel('./%s_enrich_table.xlsx'%(t), header=True, index=False)
    df.head()
    
    fig,ax=plt.subplots(figsize=(10, 15))
    x = 'log2(mean(P)/mean(N))'
    y = '-log10(pvalue)'
    df.plot(kind='scatter', x=x, y=y, ax=ax)
    ratio_max = max(df[x])
    plt.axvline(x=0, ymin=0, ymax=1, ls='--', color='grey')
    plt.axhline(y=-np.log10(0.05), xmin=0, xmax=1, ls='--', color='grey')
#     ax.set_xlim(-ratio_max-1, ratio_max+1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.title('time: %s (n=%s)'%(t, df.shape[0]))
    savefn = './%s_enrich_pvalue.pdf'%(t)
    
    df['Gene Names'] = [i.split(';')[0] for i in df['Gene names']]
    texts = []
    for x,y,t in zip(df['log2(mean(P)/mean(N))'], df['-log10(pvalue)'], df['Gene Names']):
        if y > -np.log10(0.05):
#             ax.annotate(t, (x, y), fs=3)
            texts.append(plt.text(x, y, t, fontsize=12))
#     plt.tight_layout()
    adjust_text(texts, only_move={'text': 'x'})
    
    
    plt.savefig(savefn)
    plt.close()
    
    
    
    return df, rep_N_log2_ls, rep_P_log2_ls, df[cols_keep+cols_calc]

df1, rep_N_log2_ls1, rep_P_log2_ls1, df1_save = parse_IPMASS(t='T1', mode='normal')
