import pandas as pd

# dynamic regions
dsr_egg_1cell = 'data/dynamic_region/egg_cell1/window-anno.bed'
dsr_1cell_4cell = 'data/dynamic_region/cell1_cell4/window-anno.bed'
dsr_4cell_64cell = 'data/dynamic_region/cell4_cell64/window-anno.bed'
dsr_64cell_sphere = 'data/dynamic_region/cell64_sphere/window-anno.bed'
dsr_sphere_shield = 'data/dynamic_region/sphere_shield/window-anno.bed'
way2345 = 'data/dynamic_region/way2345/way2345.bed'
way1 = 'data/dynamic_region/way1/way1.bed'


def readdsr(dsr_f):
    df = pd.read_csv(dsr_f, sep="\t", header=None, usecols=[0, 1, 2, 7])
    df.columns = ['TransID', 'start(0-based)', 'end(0-based)', 'region']
    return df


def readhotregion(f):
    df = pd.read_csv(f, sep="\t", header=None)
    df.columns = ['TransID', 'start(0-based)', 'end(0-based)']
    return df


def readAlldsr(Exwriter):
    file_ls = [dsr_egg_1cell, dsr_1cell_4cell,
               dsr_4cell_64cell, dsr_64cell_sphere, dsr_sphere_shield]
    sample_ls = ['0hpf->0.4hpf', '0.4hpf->1hpf',
                 '1hpf->2hpf', '2hpf->4hpf', '4hpf->6hpf']
    for i, j in zip(file_ls, sample_ls):
        df = readdsr(i)
        df.to_excel(Exwriter, sheet_name=j, index=False)
    return Exwriter


def dynamic_region():
    writer = pd.ExcelWriter('data/dynamic_region/all_dsr_and_hot_region.xlsx')
    writer = readAlldsr(writer)
    readhotregion(way2345).to_excel(
        writer, sheet_name='hot_dsr_way2345', index=False)
    readhotregion(way1).to_excel(
        writer, sheet_name='specific_dsr_way1', index=False)
    writer.close()


dynamic_region()


# iCLIP
h4_iCLIP = 'data/iCLIP/h4/h4_all_rep12.c6.all.bed'
h6_iCLIP = 'data/iCLIP/h6/h6_all_rep12.c6.all.bed'


def readpeaks(peak_f):
    return readdsr(peak_f)


def iCLIP():
    writer = pd.ExcelWriter('data/iCLIP/iCLIP_peaks.xlsx')
    readpeaks(h4_iCLIP).to_excel(writer, sheet_name='4hpf', index=False)
    readpeaks(h6_iCLIP).to_excel(writer, sheet_name='6hpf', index=False)
    writer.close()


iCLIP()

# motif analysis
egg_cell1_de_novo = 'data/dynamic_region/egg_cell1/utr3/egg_cell1_motif_search_site_vs_pval.txt'

def read_denovo(dn):
    df = pd.read_csv(dn, sep="\t", usecols=[2,4,5,6])
    df.columns = ['motif_seq', 'foreground_occurance', 'background_occurance', 'p_value']
    return df

def read_rbp_search(rbp):
    pass
    

##gini or mean_ic
gini_6stage = 'data/Gini.6stages.t200.null40.txt'
gini_rk33 = 'data/Gini.6stages.addRK33.t200.null40.txt'
mic_6stage = 'data/Mean_reactivity.6stages.t200.null40.txt'
mic_rk33 = 'data/Mean_reactivity.6stages.addRK33.t200.null40.txt'


def icshape():
    writer = pd.ExcelWriter('data/icSHAPE/icshape_values.6stages.xlsx')
    for s, n in zip([gini_6stage, mic_6stage], ['gini_index', 'average of reactivity']):
        df = pd.read_csv(s, index_col=0, sep="\t",)
        df = df.dropna(how='all')
        print df.shape
        df.to_excel(writer, sheet_name=n, na_rep='NULL')
    writer.close()
    writer2 = pd.ExcelWriter('data/icSHAPE/icshape_values.6stages+RK33.xlsx')
    for s, n in zip([gini_rk33, mic_rk33], ['gini_index', 'average of reactivity']):
        df = pd.read_csv(s, index_col=0, sep="\t",)
        df = df.dropna(how='all')
        print df.shape
        df.to_excel(writer2, sheet_name=n, na_rep='NULL')
    writer2.close()


# tranlation efficiency
te_rep1 = 'data/riboseq/20190415_riboseq_vs_rnaseq_control_vs_rk33_rep1.TE.txt'
te_rep2 = 'data/riboseq/20190415_riboseq_vs_rnaseq_control_vs_rk33_rep2.TE.txt'


def translation_efficiency():
    df_rep1 = pd.read_csv(te_rep1, sep="\t", index_col=0, usecols=[1, 6, 7])
    df_rep1.columns = ['control_rep1', 'RK-33_rep1']
    df_rep2 = pd.read_csv(te_rep2, sep="\t", index_col=0, usecols=[1, 6, 7])
    df_rep2.columns = ['control_rep2', 'RK-33_rep2']
    te = pd.merge(df_rep1, df_rep2, how='outer',
                  left_index=True, right_index=True)
    writer = pd.ExcelWriter('data/riboseq/translation_efficiency.xlsx')
    te.to_excel(writer, sheet_name='translation_efficiency', na_rep='/')
    writer.close()
