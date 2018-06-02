import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from sklearn.metrics import roc_curve, auc

def plot(roc_npz=None):
	if roc_npz is None:
		roc_npz = '/Share/home/zhangqf7/gongjing/zebrafish/script/RBP_prediction_wenze/HuR/human_pval0.0002_nullpct0.9_zebrafish_pval0.001_nullpct0.6_sphere/HuR_zbf_random_forest_100roc.npz'
	res = np.load(roc_npz)
	roc_auc = auc(res["arr_0"], res["arr_1"])
	plt.figure(1)
	plt.clf()
	# plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Luck', alpha=.8)
	plt.plot(res["arr_0"], res["arr_1"], lw=1, alpha=0.3, label='ROC (AUC = %0.3f)' % (roc_auc))
	#plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
	plt.xlim([0, 0.05])
	plt.ylim([0, 0.3])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic example')
	plt.legend(loc = 'best')
	#plt.show()
	#plt.draw()
	plt.savefig("%s_roc.png"%(roc_npz))

def main():
	plot()

if __name__ == '__main__':
	main()