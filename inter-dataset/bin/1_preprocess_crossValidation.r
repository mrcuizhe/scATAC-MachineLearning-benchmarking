setwd("/scATAC-machineLearning-benchmarking/inter-dataset//bin")

label_train<-read.csv('../input/pbmc_5k_atac_label_v1.txt',sep = "\t",header = TRUE)

label_test<-read.csv('../input/pbmc_5k_atac_label_nextgem.txt',sep = "\t",header = TRUE)

Cross_Validation <- function(trainLabelsPath,testLabelsPath, col_Index = 1, OutputDir){
    "
    Cross_Validation
    Function returns train and test indices for 5 folds stratified across unique cell populations,
    also filter out cell populations with less than 10 cells.
    It return a 'CV_folds.RData' file which then used as input to classifiers wrappers.
    Parameters
    ----------
    TrainLabelsPath : Cell population annotations file path for trianing (.csv).
    TestLabelsPath: Cell population annotations file path for testing (.csv)
    col_Index : column index (integer) defining which level of annotation to use,
    in case of multiple cell type annotations (default is 1)
    OutputDir : Output directory defining the path of the exported file.
    "

    trainLabels <- as.matrix(trainLabelsPath)
    trainLabels <- as.vector(trainLabels[,col_Index])
    testLabels <- as.matrix(testLabelsPath)
    testLabels <- as.vector(testLabelsPath[,col_Index])

    Removed_classes <- !(table(trainLabels) > 10)
    Cells_to_Keep_train <- !(is.element(trainLabels,names(Removed_classes)[Removed_classes]))
    trainLabels <- trainLabels[Cells_to_Keep_train]

    Removed_classes <- !(table(testLabels) > 10)
    Cells_to_Keep_test <- !(is.element(testLabels,names(Removed_classes)[Removed_classes]))
    testLabels <- testLabels[Cells_to_Keep_test]

    # Split data by cell-type
    Test_Idx<- which(testLabels != '')  #all cell
    Train_Idx<- which(trainLabels != '') #all cell
    # print(Test_Idx)
    
    #   print(Train_Idx)
    #   print(Test_Idx)
    #   print(col_Index)
    save(Train_Idx,Test_Idx,col_Index,Cells_to_Keep_train,Cells_to_Keep_test,file = paste0(OutputDir, '/CV_folds_v1Train_nextgemTest.RData'))
}

Cross_Validation(label_train,label_test, 2, "../tmp")


