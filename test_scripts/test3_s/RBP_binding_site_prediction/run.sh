nohup perl const_sample8_v2.pl HuR_motif.meme HuR_uniq_trx.txt HuR gencode.v26.transcripts.ics.std.fa 2e-3 293T.invivo.out 0.9 &
# generate training data set: positive and negative
# FIMO search potential regions based on MEME motif, if region fall into CLIP regions, then postive sample, otherwise as negative sample
# -rw-rw----+ 1 zhangqf7 zhangqf 4.3M Apr 26 10:32 HuR_pn_sample_data.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 1.8M Apr 26 10:32 HuR_scan_motif_unclip_ics.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 3.2M Apr 26 10:32 HuR_scan_motif_clip_ics.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 562K Apr 26 10:31 HuR_scan_motif_unclip.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 565K Apr 26 10:31 HuR_scan_motif_clip.txt
# -rw-rw----+ 1 zhangqf7 zhangqf  17M Apr 26 10:31 HuR_motif.meme_fimo_RBPmotif_trx.out
# -rw-rw----+ 1 zhangqf7 zhangqf 1.4M Apr 26 10:31 HuR_scan_motif.txt

nohup perl const_test7_v2.pl HuR_motif.meme HuR danRer10.refSeq.transcriptome.fa 1e-4 64cell.icshape.w200.s30.T2.t200.out 0.6 &
# generate test/prediction data set
# -rw-rw----+ 1 zhangqf7 zhangqf 1.3M Apr 26 10:36 HuR_pred_motif_data.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 1.5M Apr 26 10:36 HuR_pred_motif_ics.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 447K Apr 26 10:36 HuR_pred_motif.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 1.7M Apr 26 10:36 HuR_motif.meme_pred_RBPmotif_trx.out

nohup python random_forest_RBP2.py HuR_pn_sample_data.txt HuR_pred_motif_ics.txt HuR_zbf_random_forest 0 28 28 &
# training model, then predict 1 time
# -rw-rw----+ 1 zhangqf7 zhangqf  700 Apr 26 10:51 HuR_zbf_random_forest_important_feature.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 4.4M Apr 26 10:51 HuR_zbf_random_forest_testdata.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 882K Apr 26 10:51 HuR_zbf_random_forest_valdata.txt
# -rw-rw----+ 1 zhangqf7 zhangqf  25K Apr 26 10:51 random_forest.png
# -rw-rw----+ 1 zhangqf7 zhangqf 1.6M Apr 26 10:50 HuR_zbf_random_forest_realdata.txt

# HuR_pn_sample_data.txt 29 fields, 1-28 are features, 29 is sample label (1: postive, 0: negative)
# 28 features: p-val, icshape value of each base (7 bases), average icshape of upstream/downstream regions, variation
# [zhangqf7@loginview02 HuR]$ awk '{print NF}' HuR_pn_sample_data.txt|sort|uniq -c
#   13622 29
# HuR_pred_motif_ics.txt 36 fields, 1-8 are region info, 9-36(28) are features, predict data no sample label
# [zhangqf7@loginview02 HuR]$ awk '{print NF}' HuR_pred_motif_ics.txt|sort|uniq -c
#    4080 36

nohup python random_forest_RBP3.py HuR_pn_sample_data.txt HuR_pred_motif_ics.txt HuR_zbf_random_forest 0 28 28 &
# training model, then predict 100 times, calculate average
# -rw-rw----+ 1 zhangqf7 zhangqf  34K Apr 26 11:06 random_forest.png
# -rw-rw----+ 1 zhangqf7 zhangqf  700 Apr 26 11:05 HuR_zbf_random_forest_100important_feature.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 1.6M Apr 26 11:05 HuR_zbf_random_forest_100realdata.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 1.6K Apr 26 11:05 HuR_zbf_random_forest_100roc.npz
# -rw-rw----+ 1 zhangqf7 zhangqf 4.4M Apr 26 11:05 HuR_zbf_random_forest_100testdata.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 882K Apr 26 11:05 HuR_zbf_random_forest_100valdata.txt
# -rw-rw----+ 1 zhangqf7 zhangqf 844K Apr 26 11:05 HuR_zbf_random_forest_100vote.txt
