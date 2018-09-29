#PERL=/home/martink/bin/perl
#TABLEVIEWER_DIR=/home/martink/work/circos/svn/tools/tableviewer
TABLEVIEWER_DIR=/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer
#CIRCOS_DIR=/home/martink/work/circos/svn-tags/circos-tableviewer/
WORKING_DIR=/Share/home/zhangqf5/gongjing/software/circos-tools-0.22/tools/tableviewer2
cd $WORKING_DIR
#$PERL $TABLEVIEWER_DIR/bin/parse-table -conf etc/parse-table.conf -file uploads/table.txt -segment_order=ascii,size_desc -placement_order=row,col -interpolate_type count -col_order_row -use_col_order_row -color_source row -transparency 1 -fade_transparency 0 -ribbon_layer_order=size_asc > data/parsed.txt
$PERL $TABLEVIEWER_DIR/bin/parse-table -conf etc/parse-table.conf -file uploads/table.txt -segment_order=ascii,size_desc -placement_order=row,col -interpolate_type count -col_order_row -use_col_order_row -col_color_row -use_col_color_row -ribbon_layer_order=size_asc > data/parsed.txt
cat data/parsed.txt | $TABLEVIEWER_DIR/bin/make-conf -dir data
circos -param random_string=zgvickusamp -conf etc/circos.conf
cd -

#cd $WORKING_DIR; circos -param random_string=zgvickusamp -conf ${WORKING_DIR}/etc/circos.conf 2>&1 > ${WORKING_DIR}/results/report.txt
#tar cvfz circos-tableviewer-zgvickusamp.tar.gz *
