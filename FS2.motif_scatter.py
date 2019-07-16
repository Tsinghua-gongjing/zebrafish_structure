def plot_motif_logo(motif_ls, site_num_ls, pvalue_ls, tomtom_file, rbp_name):
    dd = pd.DataFrame({'site_num':site_num_ls, 'pvalue':pvalue_ls})
    dd['-log10(pvalue)'] = -np.log10(dd['pvalue'])
    dd['log2(site_num)'] = np.log2(dd['site_num'])
    dd['motif'] = motif_ls
    print dd
    sns.set_style('ticks')
    fig,ax = plt.subplots(figsize=(25,20))
    dd.plot(kind='scatter', x='log2(site_num)', y='-log10(pvalue)', ax=ax)
    ax.set_xlim(0,)
    
    tomtom = read_tomtom(tomtom_file, rbp_name)
    
    texts = []
    for x,y,t in zip(dd['log2(site_num)'], dd['-log10(pvalue)'], dd['motif']):
        if tomtom.has_key(t):
            RBP = t + '\n'+'\n'.join(tomtom[t])
        else:
            RBP = t
        if t in ['TDTWV', 'TDWDT', 'AMACAMAC', 'TAANC', 'ACAYT']:
            texts.append(plt.text(x, y, RBP, fontsize=12))
    adjust_text(texts, only_move={'text': 'xy'})
    
    plt.title(tomtom_file.split('/')[-2])
    
    plt.tight_layout()
    savefn = './%s/%s.pdf'%(tomtom_file.split('/')[-3], tomtom_file.split('/')[-2])
    plt.savefig(savefn)
    
motif_df = read_dmeme('./d10/way2345/dreme.txt')
motif_ls = list(motif_df[2])
site_num_ls2 = map(int, list(motif_df[3]))
pvalue_ls2 = [1.0e-100 if i < 1.0e-100 else i for i in map(float, list(motif_df[5]))]
tomtom_file = './d10/way2345/tomtom.tsv'

plot_motif_logo(motif_ls, site_num_ls2, pvalue_ls2, tomtom_file, rbp_name)