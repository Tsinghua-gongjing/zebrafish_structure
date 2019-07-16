#!/usr/bin/python

import re
import sys

import matplotlib as mpl
mpl.use('Agg')
import time

import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import make_blobs
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import log_loss
from random import shuffle

from itertools import cycle
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp



if __name__ == '__main__':
    testdata = sys.argv[1]
    realdata = sys.argv[2]
    relaspe = sys.argv[3]
    stacord = sys.argv[4]
    endcord = sys.argv[5]
    labelcord = sys.argv[6]

    np.random.seed(0)

    X = []
    y = []
    strline = []
    infile = open(testdata,'rb')
    for line in infile:
        strline.append(line)

    infile.close()
    
    Real_ind = []
    rdata = []
    vote_re = []
    final_re = []
    infile = open(realdata,'rb')
    for line in infile:
        line = line.strip('\n')
        sent = line.split('\t')
        Real_ind.append(sent)
        vote_re.append(sent[0])
        final_re.append(0)
        xi = sent[8 + int(stacord): 8 + int(endcord)]
        xi = np.asarray(xi,dtype = float)
        rdata.append(xi)
        
    Real_ind = np.asarray(Real_ind)
    rdata = np.asarray(rdata)
    vote_re = np.asarray(vote_re)
    final_re = np.asarray(final_re)
    
    vote_num = 0
    total_auc = 0.0
    roc_auc_max = 0.0
    while vote_num < 100:
        shuffle(strline)
        Test_ind = []
        X = []
        y = []
        datanum = 0
        for line in strline:
            datanum = datanum + 1
            line = line.strip('\n')
            sent = line.split('\t')
            Test_ind.append(sent)
            xi = sent[int(stacord):int(endcord)]
            yi = int(sent[int(labelcord)])
            xi = np.asarray(xi,dtype = float)
            X.append(xi)
            y.append(yi)

        X = np.asarray(X)
        y = np.asarray(y)
        Test_ind = np.asarray(Test_ind)

        X_train, y_train = X[:int(datanum*0.6)], y[:int(datanum*0.6)]
        X_valid, y_valid = X[int(datanum*0.6):int(datanum*0.8)], y[int(datanum*0.6):int(datanum*0.8)]
        X_train_valid, y_train_valid = X[:int(datanum*0.8)], y[:int(datanum*0.8)]
        X_test, y_test = X[int(datanum*0.8):], y[int(datanum*0.8):]

        clf = RandomForestClassifier(n_estimators=100)
        clf.fit(X_train_valid, y_train_valid)
        clf_probs = clf.predict_proba(X_test)
        score = log_loss(y_test, clf_probs)
        print clf.score(X_test, y_test)
        x_pro = clf.predict_proba(X)
        Test_res = np.c_[Test_ind, x_pro]
        val_res = np.c_[X_test, y_test]
        val_res = np.c_[val_res, clf_probs]

        fpr, tpr, thresholds = roc_curve(y_test, clf_probs[:, 1])
        roc_auc = auc(fpr, tpr)
        print roc_auc
        if roc_auc > 0.65:
            if roc_auc > total_auc/(vote_num + 1):
                roc_auc_max = roc_auc
                fpr_max = fpr
                tpr_max = tpr
            total_auc = total_auc + roc_auc
            vote_num = vote_num + 1
            yy_pre = clf.predict(rdata)
            yy_pro = clf.predict_proba(rdata)
            final_re = final_re + yy_pre
            vote_re = np.c_[vote_re,yy_pre]
            
        
    final_re = final_re.astype('float')
    final_re = final_re/vote_num
    rdata_res = np.c_[Real_ind,yy_pre]
    rdata_res = np.c_[rdata_res,yy_pro]
    rdata_res = np.c_[rdata_res,final_re]
    print total_auc/vote_num
    rdata_file = relaspe + "_100realdata.txt"
    tdata_file = relaspe + "_100testdata.txt"
    valdata_file = relaspe + "_100valdata.txt"
    feature_file = relaspe + "_100important_feature.txt"
    vote_file = relaspe + "_100vote.txt"
    roc_file = relaspe + "_100roc.npz"
    np.savez(roc_file, fpr_max, tpr_max)
    np.savetxt(vote_file, vote_re, fmt='%s', newline='\n')
    np.savetxt(rdata_file, rdata_res, fmt='%s', newline='\n')
    np.savetxt(tdata_file, Test_res, fmt='%s', newline='\n')
    np.savetxt(valdata_file, val_res, fmt='%s', newline='\n')
    np.savetxt(feature_file, clf.feature_importances_)

    # Compute ROC curve and ROC area for each class
    
    res = np.load(roc_file)
    roc_auc = auc(res["arr_0"], res["arr_1"])
    plt.figure(1)
    plt.clf()
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Luck', alpha=.8)
    plt.plot(res["arr_0"], res["arr_1"], lw=1, alpha=0.3, label='ROC (AUC = %0.3f)' % (roc_auc))
    #plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc = 'best')
    #plt.show()
    #plt.draw()
    plt.savefig("random_forest.png")



