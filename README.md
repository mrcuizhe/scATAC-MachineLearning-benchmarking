# scATAC-MachineLearning-benchmarking

Single-cell assay for transposase accessible chromatin using sequencing(scATAC-seq) is rapidly advancing our understanding of the cellular composition of complex tissues and organisms. The similarity of data structure and feature between scRNA-seq and scATAC-seq makes it feasible to identify the cell types in scATAC-seq through traditional supervised machine learning methods, which is widely proved reliable and used as underlying method in scRNA-seq classification analysis. Here, we evaluated 6 popular machine learning methods for classification in scATAC-seq. The performance of the methods is evaluated using 4 publicly available single cell ATAC-seq datasets of different tissues and sizes and technologies. We evaluated these methods using intra-datasets experiments of 5-folds cross validation based on accuracy, percentage of correctly predicted cells. We found that these methods may perform well in some types of cells in single dataset, but the overall results are not as well as in scRNA-seq analysis. The SVM method has overall the best performance across all experiments.

## Description
We provide all the scripts to run and evaluate all machine laerning methods, and to reproduce the results introduced in the paper.

All processed matrix used in paper are also provided.

1. Each folder stored the codes, input and output for a experiment
2. The bin folder stores the benchmarking and evaluation code
3. The input folder stores the bin-cell matrix and labels for dataset
4. The tmp folder stores the temporary file of experiment, especially the grouping information for 5-fold cross validation
5. The output file stores the results of each experiment and its evaluation results


## Citation
XXX

## Credits
Zhe Cui, Liran Juan, Tao Jiang, Bo Liu, Tianyi Zang and Yadong Wang

## Email
 cuizhe@hit.edu.cn, mrcuizhe@gmail.com