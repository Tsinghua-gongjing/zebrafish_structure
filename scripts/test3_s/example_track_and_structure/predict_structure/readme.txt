运行：
python structure_predict.py icshape.out transcript_id start_coordinate end_coordinate reference_fasta

例子：
python structure_predict.py /Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/1cell.icshape.w200.s30.T2.t200.out NM_001123007 100 130 /Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.fa

输出：
predict structure: GGGGTTGTGAACTCAATGAATCTAATAAAG .(((((...)))))................

note:
coordinate: 0-based
rnastructure: 软件需要安装，使用了其中的ct2dot, Fold-smp命令

[zhangqf7@bnode02 script]$ which ct2dot
/Share/home/zhangqf/usr/rnastructure/ct2dot
[zhangqf7@bnode02 script]$ which Fold-smp
/Share/home/zhangqf/usr/rnastructure/Fold-smp
