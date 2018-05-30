import subprocess, os
import sys

def tableview(txt=None, parse_table=None, make_conf=None, conf_dir=None, circos_conf=None, save_dir=None, parsed_conf=None):
	if parse_table is None:
		parse_table = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer/bin/parse-table'
	if make_conf is None:
		make_conf = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer/bin/make-conf'
	if conf_dir is None:
		conf_dir = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer2/data'
	if circos_conf is None:
		circos_conf = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer2/etc/circos.conf'
	if save_dir is None:
		save_dir = os.path.dirname(txt)
	if parsed_conf is None:
		#parsed_conf = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer/samples/parse-table-02a.conf'
		parsed_conf = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer2/etc/parse-table.conf'
	print "[tableview start] file: %s"%(txt)
	subprocess.call(["cat {txt} | {parse_table} -conf {parsed_conf} -segment_order=ascii,size_desc -placement_order=row,col -interpolate_type count -col_order_row -use_col_order_row -col_color_row -use_col_color_row -ribbon_layer_order=size_asc | {make_conf} -dir {conf_dir}".format(txt=txt, parse_table=parse_table, parsed_conf=parsed_conf, make_conf=make_conf, conf_dir=conf_dir)],shell=True)
	#subprocess.call(["cat {txt} | {parse_table} -conf {parsed_conf} | {make_conf} -dir {conf_dir}".format(txt=txt, parse_table=parse_table, parsed_conf=parsed_conf, make_conf=make_conf, conf_dir=conf_dir)],shell=True)
	file_png = txt.split('/')[-1].replace('txt', 'png')
	subprocess.call(["circos -conf {circos_conf} -outputdir {save_dir} -outputfile {file_png} -param random_string=zgvickusamp| grep created".format(circos_conf=circos_conf, save_dir=save_dir, file_png=file_png)],shell=True)
	print "[tableview end] file: %s"%(txt)

def main():
	#tableview(txt='/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer/samples/RRI_union_deduplicates.split_full.type_dis.txt')
        if len(sys.argv) == 1:
            txt = '/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer2/samples/RRI_union_deduplicates.split_full.type_dis.txt'
        else:
            txt = sys.argv[1]
        tableview(txt)

if __name__ == '__main__':
	main()
