import pandas as pd
from scipy import stats

sns.set(style="ticks")
sns.set_context("poster")

gini = './RPKM_combine.merge.t200.gini.null40.txt'
TE_bartel = './TE.matrix' 

df_gini = pd.read_csv(gini, header=0, index_col=0, sep='\t')
print df_gini.head()

df_TE_dev = pd.read_csv(TE_bartel, header=0, sep='\t')
df_TE_dev['id'] = df_TE_dev.index
print df_TE_dev.head() 

def plot_gini_TE_ratio(TE_time_col1='256cell', TE_time_col2='dome', stage_col1='64cell(transcript)', stage_col2='sphere(transcript)'):
    df_gini['id'] = df_gini.index

    df_gini_2 = df_gini[['id', stage_col1, stage_col2]]
    df_gini_2 = df_gini_2[(df_gini_2[stage_col1]>0) & (df_gini_2[stage_col2]>0)]
    df_gini_2['%s/%s'%(stage_col2, stage_col1)] = df_gini_2[stage_col2]/df_gini_2[stage_col1]
#     print df_gini_2.head()

    df_TE_dev_2 = df_TE_dev[['id', TE_time_col1, TE_time_col2]]
    df_TE_dev_2 = df_TE_dev_2[(df_TE_dev_2[TE_time_col1]>0) & (df_TE_dev_2[TE_time_col2]>0)]
    df_TE_dev_2['%s/%s'%(TE_time_col2, TE_time_col1)] = df_TE_dev_2[TE_time_col2] / df_TE_dev_2[TE_time_col1]
    df_TE_dev_2['%s-%s'%(TE_time_col2, TE_time_col1)] = df_TE_dev_2[TE_time_col2] - df_TE_dev_2[TE_time_col1]

    df_gini_TE = pd.merge(df_gini_2, df_TE_dev_2, on='id')
#     print df_gini_TE.head(), df_gini_TE.describe()
    print df_gini_TE.shape

#     sns.jointplot(y='%s-%s'%(TE_time_col2, TE_time_col1), x='%s/%s'%(stage_col2, stage_col1), 
#                   data=df_gini_TE, stat_func=stats.pearsonr, size=10,)
    sns.jointplot(y='%s-%s'%(TE_time_col2, TE_time_col1), x='%s/%s'%(stage_col2, stage_col1), 
                  data=df_gini_TE, kind="kde")
    plt.savefig('./png/TE.%s_%s.pdf'%(stage_col1, stage_col2))
    plt.close()
    
    
# plot_gini_TE_ratio(TE_time_col1='256cell', TE_time_col2='dome', stage_col1='64cell(transcript)', stage_col2='sphere(transcript)')
# plot_gini_TE_ratio(TE_time_col1='256cell', TE_time_col2='dome', stage_col1='64cell(UTR5)', stage_col2='sphere(UTR5)')
# plot_gini_TE_ratio(TE_time_col1='256cell', TE_time_col2='dome', stage_col1='64cell(CDS)', stage_col2='sphere(CDS)')
# plot_gini_TE_ratio(TE_time_col1='256cell', TE_time_col2='dome', stage_col1='64cell(UTR3)', stage_col2='sphere(UTR3)')


# plot_gini_TE_ratio(TE_time_col1='dome', TE_time_col2='shield', stage_col1='sphere(transcript)', stage_col2='shield(transcript)')
# plot_gini_TE_ratio(TE_time_col1='dome', TE_time_col2='shield', stage_col1='sphere(UTR5)', stage_col2='shield(UTR5)')
# plot_gini_TE_ratio(TE_time_col1='dome', TE_time_col2='shield', stage_col1='sphere(CDS)', stage_col2='shield(CDS)')
# plot_gini_TE_ratio(TE_time_col1='dome', TE_time_col2='shield', stage_col1='sphere(UTR3)', stage_col2='shield(UTR3)')


plot_gini_TE_ratio(TE_time_col1='hour2', TE_time_col2='hour4', stage_col1='64cell(transcript)', stage_col2='sphere(transcript)')
plot_gini_TE_ratio(TE_time_col1='hour2', TE_time_col2='hour4', stage_col1='64cell(UTR5)', stage_col2='sphere(UTR5)')
plot_gini_TE_ratio(TE_time_col1='hour2', TE_time_col2='hour4', stage_col1='64cell(CDS)', stage_col2='sphere(CDS)')
plot_gini_TE_ratio(TE_time_col1='hour2', TE_time_col2='hour4', stage_col1='64cell(UTR3)', stage_col2='sphere(UTR3)')


plot_gini_TE_ratio(TE_time_col1='hour4', TE_time_col2='hour6', stage_col1='sphere(transcript)', stage_col2='shield(transcript)')
plot_gini_TE_ratio(TE_time_col1='hour4', TE_time_col2='hour6', stage_col1='sphere(UTR5)', stage_col2='shield(UTR5)')
plot_gini_TE_ratio(TE_time_col1='hour4', TE_time_col2='hour6', stage_col1='sphere(CDS)', stage_col2='shield(CDS)')
plot_gini_TE_ratio(TE_time_col1='hour4', TE_time_col2='hour6', stage_col1='sphere(UTR3)', stage_col2='shield(UTR3)')

plot_gini_TE_ratio(TE_time_col1='hour2', TE_time_col2='hour6', stage_col1='64cell(transcript)', stage_col2='shield(transcript)')
plot_gini_TE_ratio(TE_time_col1='hour2', TE_time_col2='hour6', stage_col1='64cell(UTR5)', stage_col2='shield(UTR5)')
plot_gini_TE_ratio(TE_time_col1='hour2', TE_time_col2='hour6', stage_col1='64cell(CDS)', stage_col2='shield(CDS)')
plot_gini_TE_ratio(TE_time_col1='hour2', TE_time_col2='hour6', stage_col1='64cell(UTR3)', stage_col2='shield(UTR3)')
