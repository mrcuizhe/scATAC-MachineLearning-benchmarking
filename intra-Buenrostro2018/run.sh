#!/bin/bash

cd /scATAC-machineLearning-benchmarking/intra-Buenrostro2018/bin

Rscript 1_preprocess_crossValidation.r

python3 2_run_RF.py
python3 2_run_NMC.py
python3 2_run_LDA.py
python3 2_run_SVM.py
python3 2_run_DT.py
python3 2_run_KNN50.py
python3 2_run_KNN9.py

Rscript 3_evaluate_RF.r 
Rscript 3_evaluate_NMC.r
Rscript 3_evaluate_LDA.r 
Rscript 3_evaluate_KNN50.r 
Rscript 3_evaluate_KNN9.r 
Rscript 3_evaluate_DT.r 
Rscript 3_evaluate_SVM.r 