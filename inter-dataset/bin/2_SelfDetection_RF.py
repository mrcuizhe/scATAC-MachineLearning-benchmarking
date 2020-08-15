#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
from sys import argv
from pathlib import Path
import numpy as np
import pandas as pd
import time as tm
from sklearn.ensemble import RandomForestClassifier
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter
from scipy import sparse


# In[3]:


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


# In[6]:


def run_RF(trainDataPath,testDataPath, trainLabelsPath,testLabelsPath, CV_RDataPath, OutputDir):
    '''
    run baseline classifier: RF
    Wrapper script to run an RF classifier on a benchmark dataset with 5-fold cross validation,
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

    tokeep_train = np.array(ro.r['Cells_to_Keep_train'], dtype = 'bool')
    tokeep_test = np.array(ro.r['Cells_to_Keep_test'], dtype = 'bool')
    col = np.array(ro.r['col_Index'], dtype = 'int')
    col = col - 1
    test_ind =ro.r['Test_Idx']
    train_ind =ro.r['Train_Idx']

    # read the data
    data_train=ro.r['readRDS'](trainDataPath)
    data_train=pd.DataFrame(dgc_to_csr(data_train).toarray(),dtype='uint8').T
    
    data_test=ro.r['readRDS'](testDataPath)
    data_test=pd.DataFrame(dgc_to_csr(data_test).toarray(),dtype='uint8').T
    
    labels_train = pd.read_csv(trainLabelsPath, header=0,index_col=None, sep='\t', usecols = col)
    labels_test = pd.read_csv(testLabelsPath, header=0,index_col=None, sep='\t', usecols = col)
    
#     print(len(data))
#     print(labels)
#     print(len(tokeep))
    labels_train = labels_train.iloc[tokeep_train]
    labels_test = labels_test.iloc[tokeep_test]
    data_train = data_train.iloc[tokeep_train]
    data_test = data_test.iloc[tokeep_test]

    # normalize data
    data_train = np.log1p(data_train)
    data_test = np.log1p(data_test)

    Classifier = RandomForestClassifier(n_estimators = 50)

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []
    
    test_ind = np.array(test_ind, dtype = 'int') - 1
    train_ind = np.array(train_ind, dtype = 'int') - 1

    train=data_train.iloc[train_ind]
    test=data_test.iloc[test_ind]
    y_train=labels_train.iloc[train_ind]
    y_test=labels_test.iloc[test_ind]

    print(test_ind)
    start=tm.time()
    Classifier.fit(train, y_train)
    tr_time.append(tm.time()-start)

    start=tm.time()
    predicted = Classifier.predict(test)
    ts_time.append(tm.time()-start)

    truelab.extend(y_test.values)
    pred.extend(predicted)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)
#     print(len(tr_time))

    OutputDir = Path(OutputDir)
    os.makedirs(Path(OutputDir),exist_ok=True)
    truelab.to_csv(str(OutputDir / Path("RF_true.csv")),
                   index = False)
    pred.to_csv(str(OutputDir / Path("RF_pred.csv")),
                index = False)
    tr_time.to_csv(str(OutputDir / Path("RF_training_time.csv")),
                   index = False)
    ts_time.to_csv(str(OutputDir / Path("RF_test_time.csv")),
                   index = False)

run_RF("../input/10xpbmc5k-snap-full_v1Train.rds","../input/10xpbmc5k-snap-full_nextgemTest.rds", "../input/pbmc_5k_atac_label_v1.txt","../input/pbmc_5k_atac_label_nextgem.txt", "../tmp/CV_folds_v1Train_nextgemTest.RData", "../output")





