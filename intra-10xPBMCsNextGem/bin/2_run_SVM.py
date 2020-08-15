#!/usr/bin/env python
# coding: utf-8

import os
from sys import argv
from pathlib import Path
import numpy as np
import pandas as pd
import time as tm
from sklearn.svm import LinearSVC
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter
from scipy import sparse

ro.r("library(Matrix)")

def dgc_to_csr(r_dgc):
    """Convert (and transpose) a dgCMatrix from R to a csr_matrix in python
    """
    with localconverter(ro.default_converter + pandas2ri.converter):
        X = sparse.csr_matrix(
                (
                    r_dgc.slots["x"], 
                    r_dgc.slots["i"], 
                    r_dgc.slots["p"]
                ),
                shape=tuple(ro.r("dim")(r_dgc))[::-1]
            )
    return X

def run_SVM(DataPath, LabelsPath, CV_RDataPath, OutputDir):
    '''
    run baseline classifier: SVM
    Wrapper script to run an SVM classifier with a linear kernel on a benchmark dataset with 5-fold cross validation,
    outputs lists of true and predicted cell labels as csv files, as well as computation time.
    Parameters
    ----------
    DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes
    as row names and gene names as column names.
    LabelsPath : Cell population annotations file path (.csv).
    CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
    OutputDir : Output directory defining the path of the exported file.
    '''

    # read the Rdata file
    ro.r['load'](CV_RDataPath)

    nfolds = np.array(ro.r['n_folds'], dtype = 'int')
    tokeep = np.array(ro.r['Cells_to_Keep'], dtype = 'bool')
    col = np.array(ro.r['col_Index'], dtype = 'int')
    col = col - 1
    test_ind =ro.r['Test_Idx']
    train_ind =ro.r['Train_Idx']

    # read the data
    data=ro.r['readRDS'](DataPath)
    data=pd.DataFrame(dgc_to_csr(data).toarray()).T
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep='\t', usecols = col)
    
#     print(len(data))
#     print(labels)
#     print(len(tokeep))
    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]

    # normalize data
    data = np.log1p(data)

    Classifier = LinearSVC()

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []
    
    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        start=tm.time()
        Classifier.fit(train, y_train)
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = Classifier.predict(test)
        ts_time.append(tm.time()-start)

        truelab.extend(y_test.values)
        pred.extend(predicted)
        print(len(pred))

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)
#     print(len(tr_time))

    OutputDir = Path(OutputDir)
    os.makedirs(Path(OutputDir),exist_ok=True)
    truelab.to_csv(str(OutputDir / Path("SVM_true.csv")),
                   index = False)
    pred.to_csv(str(OutputDir / Path("SVM_pred.csv")),
                index = False)
    tr_time.to_csv(str(OutputDir / Path("SVM_training_time.csv")),
                   index = False)
    ts_time.to_csv(str(OutputDir / Path("SVM_test_time.csv")),
                   index = False)

run_SVM("../output/10xpbmc5k-snap-full.rds", "../input/pbmc_5k_atac_label_nextgem.txt", "../tmp/CV_folds.RData", "../output/")



