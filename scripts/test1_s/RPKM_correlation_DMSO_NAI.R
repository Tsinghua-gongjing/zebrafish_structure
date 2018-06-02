df_to_corr_matrix = function(df, feature) {
	df1 = acast(df, sample1 ~ sample2, value.var=feature, fill=0)
	df2 = acast(df, sample2 ~ sample1, value.var=feature, fill=0)
	df3 = df1 + df2 + diag(x=-1, 32, 32)
	return(df3)
}

merge_plot = function(plotData, title){
	pheatmap(plotData, lwd=2, color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")), bias=2)(100), clustering_method="complete", main=title, display_numbers=F, number_format="%.2f", fontsize_number = 0.4 * 10, cluster_rows = FALSE, cluster_cols = FALSE)
}

setwd("/Users/soul/baiduyun/research/zhanglab/1-YANG_lab_cooperation_projects/zebrash/data_analysis/paper_figures/icSHAPE/RPKM/")

library(reshape2)
library(pheatmap)
library(RColorBrewer)
f = 'stat.spearman.txt'

df = read.table(f, header=T, sep='\t')

print(head(df))

feature_ls = c('pearson')

#RT_files = '/Share/home/zhangqf7/gongjing/zebrafish/data/RT'
#cols_all = list.files(RT_files)

#cols_DMSO = list.files(RT_files, pattern='DMSO*')
#cols_NAI = list.files(RT_files, pattern='NAI*')

cols_sort = c('DMSO_egg_rep1', 'DMSO_egg_rep2', 'DMSO_1cell_rep1', 'DMSO_1cell_rep2', 'DMSO_4cell_rep2', 'DMSO_4cell_rep1', 'DMSO_64cell_rep1', 'DMSO_64cell_rep2', 
	          'DMSO_sphere_rep1', 'DMSO_sphere_rep2', 'DMSO_shield_rep1', 'DMSO_shield_rep2',
	          'NAI_egg_rep1', 'NAI_egg_rep2', 'NAI_1cell_rep1', 'NAI_1cell_rep2', 'NAI_4cell_rep2', 'NAI_4cell_rep1', 'NAI_64cell_rep1', 'NAI_64cell_rep2', 
	           'NAI_sphere_rep1', 'NAI_sphere_rep2', 'NAI_shield_rep1', 'NAI_shield_rep2')

cols_dmso = c('DMSO_egg_rep1', 'DMSO_egg_rep2', 'DMSO_1cell_rep1', 'DMSO_1cell_rep2', 'DMSO_4cell_rep1', 'DMSO_4cell_rep2', 'DMSO_64cell_rep1', 'DMSO_64cell_rep2', 
               'DMSO_sphere_rep1', 'DMSO_sphere_rep2', 'DMSO_shield_rep1', 'DMSO_shield_rep2')
cols_nai = c('NAI_egg_rep1', 'NAI_egg_rep2', 'NAI_1cell_rep1', 'NAI_1cell_rep2', 'NAI_4cell_rep1', 'NAI_4cell_rep2', 'NAI_64cell_rep1', 'NAI_64cell_rep2', 
              'NAI_sphere_rep1', 'NAI_sphere_rep2', 'NAI_shield_rep1', 'NAI_shield_rep2')

for(feature in feature_ls){

		print(paste('plot: ', feature, sep=' '))

		df_plot_corr_marix = df_to_corr_matrix(df, feature)
		print(df_plot_corr_marix)

#		merge_plot(plotData=df_plot_corr_marix[cols_all, cols_all], title=paste('RPKM', 'all', sep=' '))
#		merge_plot(plotData=df_plot_corr_marix[cols_NAI, cols_NAI], title=paste('RPKM', 'NAI', sep=' '))
#		merge_plot(plotData=df_plot_corr_marix[cols_DMSO, cols_DMSO], title=paste('RPKM', 'DMSO', sep=' '))
                df_to_plot = df_plot_corr_marix[cols_sort, cols_sort]
                cols_names = c('DMSO_egg_rep1', 'DMSO_egg_rep2', 'DMSO_1cell_rep1', 'DMSO_1cell_rep2', 'DMSO_4cell_rep1', 'DMSO_4cell_rep2', 'DMSO_64cell_rep1', 'DMSO_64cell_rep2', 
	           'DMSO_sphere_rep1', 'DMSO_sphere_rep2', 'DMSO_shield_rep1', 'DMSO_shield_rep2', 
	          'NAI_egg_rep1', 'NAI_egg_rep2', 'NAI_1cell_rep1', 'NAI_1cell_rep2', 'NAI_4cell_rep1', 'NAI_4cell_rep2', 'NAI_64cell_rep1', 'NAI_64cell_rep2', 
	           'NAI_sphere_rep1', 'NAI_sphere_rep2', 'NAI_shield_rep1', 'NAI_shield_rep2')
                colnames(df_to_plot) = cols_names
                rownames(df_to_plot) = cols_names
		#merge_plot(plotData=df_to_plot, title=paste('RPKM Correlation', 'all (sort by time)', sep=' '))
		df_dmso_plot = df_to_plot[cols_dmso, cols_dmso]
		colnames(df_dmso_plot) = c('0h.p.f rep1', '0h.p.f rep2', '0.2h.p.f rep1', '0.2h.p.f rep2', '1h.p.f rep1', '1h.p.f rep2', '2h.p.f rep1', 
		                          '2h.p.f rep2', '4h.p.f rep1', '4h.p.f rep2',  '6h.p.f rep1', '6h.p.f rep2')
		rownames(df_dmso_plot) = colnames(df_dmso_plot)
		df_nai_plot = df_to_plot[cols_nai, cols_nai]
		colnames(df_nai_plot) = colnames(df_dmso_plot)
		rownames(df_nai_plot) = colnames(df_dmso_plot)
		pdf(paste('RPKMPairWiseCorrelationPlot_DMSO.', format(Sys.time(), "%Y%m%dhr%Hm%M",), ".pdf"), onefile = FALSE)
		merge_plot(plotData=df_dmso_plot, title = paste('RPKM Correlation', 'DMSO', sep=' '))
		dev.off()
		pdf(paste('RPKMPairWiseCorrelationPlot_NAI.', format(Sys.time(), "%Y%m%dhr%Hm%M",), ".pdf"), onefile = FALSE)
		merge_plot(plotData=df_nai_plot, title = paste('RPKM Correlation', 'NAI', sep=' '))
		dev.off()
}




