# DeepCirCode 
For circular RNA (circRNA) back-splicing prediction, DeepCirCode utilizes a convolutional neural network (CNN) with one-hot encoded nucleotide sequence as the input. Relevant features can be automatically extracted by DeepCirCode. Each kernel of the first CNN layer of DeepCirCode can be regarded as a motif scanner. In order to make these DeepCirCode-learnt features biologically interpretable, we have also implemented a visualization method to represent these features as sequence motifs. The DeepCirCode method has been applied to three different species, including human, mouse and fruit fly, for circRNA back-splicing prediction. Comparision of DeepCirCode-learnt motifs from these species has also been performed.  

This documentation is part of the supplementary information release for DeepCirCode. For details of this work, users can refer to our paper "**Deep Learning of the Back-splicing Code for Circular RNA Formation**" (J. Wang  and L. Wang, 2019). 
# Requirements 
DeepCirCode is an R package with one necessary function written in Python. To use DeepCirCode package, R >= 3.5.3 is required. The versions of other packages are suggestive ones. 

- R 3.5.3 + 64-bit; 

- Python 2.7 + 64-bit; 

- Keras 2.2.4 in R and Python; 

- Tensorflow 1.10 in R and Python; 

- ROCR 1.0-7 in R

- reticulate 1.11.1 in R

- numpy in Python; 
# Installing DeepCirCode 
``` 
# Download "DeepCirCode" package using devtools from github, and load DeepCirCode: 

devtools::install_github("BioDataLearning/DeepCirCode") 

library("DeepCirCode") 

# users can check whether DeepCirCode has been loaded by using: 

devtools::loaded_packages() 
``` 
# Getting the training and test sets 
The training and test sets for all the three species (human, mouse and fruit fly) have been included in the "DeepCirCode" package. They are named as "HumanTrain", "HumanTest", "MouseTrain", "MouseTest", "FlyTrain" and "FlyTest". Users can load a dataset by using the function **data(**dataset_Name**)**, such as: 
``` 
data(HumanTrain)  # datasets are lazy loaded as variable “HumanTrain” until being used
data(HumanTest) 
``` 
All the datasets are in the form of c(instance_number,3), with each of the three variables representing "**label**" (0 for negative, 1 for positive),"**RNA_input_seq**" (the raw 200 nucleotide RNA sequences, e.g. AUCC...GACG) and "**encoded_seq**" (one-hot encoded RNA sequence of 800 characters in length, e.g. 100000010001...00010010), which can be checked using: 
``` 
dim(HumanTrain) 
colnames(HumanTrain) 

HumanTrain$label[1] 
HumanTrain$RNA_input_seq[1] 
HumanTrain$encoded_seq[1] 
``` 
# Converting one-hot encoded sequence into 3D matrix input for DeepCirCode 
For convenience, all the datasets of DeepCirCode have been included in this package. This section provides users some help for using DeepCirCode package to prepare their own datasets. Users may also jump to next section for dataset loading and model construction. 

The first layer of DeepCirCode is a 1D convolution layer, which requires the input_shape = c(200,4) for sequences of 200-nt with 4 channels (A,T/U,G,C). **x_train** (features of training set) and **x_test** (features of test set) should be converted into matrix of the shape c(instance_number, 200,4); **y_train** (label of training set) and **y_test** (label of test set) should be converted into matrix of the shape c(instance_number, 2) with **y[ ,1]** represents the negative class and **y[ ,2]** represents the positive class. 

To convert one-hot encoded string to 3D matrix, use DeepCirCode function **convStringToMatrix()**. Here, use human dataset as an example: 
``` 
x_train <- convStringToMatrix(HumanTrain$encoded_seq) 
x_test <- convStringToMatrix(HumanTest$encoded_seq) 

# This process may take several minutes depends on the size of dataset 
# x_train[1,2,] will print out the one-hot encoded vector of the 2nd nucleotide of HumanTrain$RNA_input_seq[1] 
``` 
To convert label(1/0) into categorical matrix, use Keras function **to_categorical()**: 
``` 
library(keras) 
library(tensorflow) 

y_train <- to_categorical(HumanTrain$label,2) 
y_test <- to_categorical(HumanTest$label,2) 
``` 
# Testing DeepCirCode with test sets 
**Firstly**, load required packages 
``` 
# Loading packages: 

library(keras) 
library(tensorflow) 
library(ROCR) 
``` 
**Secondly**, load test sets in matrix: 

All the x_train, x_test, y_train, y_test datasets have been packaged in DeepCirCode. Instead of re-creating datasets by following the above mentioned steps, users can directly load these input datasets by function **data(** dataset_Name **)**: 
``` 
# Loading datasets, using human datasets as an example: 

data(x_test_human) # data has been loaded and stored in the variable "x_test_human" 
data(y_test_human) 
``` 
**Thirdly**, applying DeepCirCode model on the test sets: 

The DeepCirCode models for all three species have been packaged under the main directory of DeepCirCode named "**DeepCirCode_bestmodel_human.hdf5**", "**DeepCirCode_bestmodel_mouse.hdf5**" and "**DeepCirCode_bestmodel_fly.hdf5**". To perform the prediction of DeepCirCode on the test sets, using DeepCirCode function **DeepCirCode_test()** which takes arguments the x_test, y_test, and species(="human" or "mouse" or "fly") such as: 
``` 
DeepCirCode_test(x_test_human, y_test_human, "human") # don't forget the "" for species argument 

# The DeepCirCode model architecture, prediction confusion matrix, loss, accuracy, ROC curve and ROC_AUC value will be printed. 
``` 
# Motif visualization 
We applied the strategy as described by [DanQ](https://github.com/uci-cbcl/DanQ/issues/9) for motif visualization written in python script, named "**getMotifs.py**" under the directory of "DeepCircode". "getMotifs.py" defines a function **DeepCirCode_getMotifs()** which takes input x_train, y_train, x_test, y_test, and outputs a "**DeepCirCode_Position_Frequency_Matrix.txt**" file containing the position frequency matrix for each kernel, which can be directly uploaded to [TOMTOM](http://meme-suite.org/tools/tomtom) for motif comparision and visualization. 

Here, we applied R package "reticulate" to run this python script in R interface. 

**Firstly**, load package and "getMotifs.py": 
``` 
library(reticulate) 
path <- system.file("getMotifs.py", package = "DeepCirCode")
source_python(path) 
``` 
**Secondly**, run DeepCirCode_getMotifs(): 

Since motif visualization needs the information of the first layer of DeepCirCode which is not stored in a saved model, we need to train DeepCirCode freshly. So, DeepCirCode_getMotifs() will train DeepCirCode and then visualize the kernels of the first CNN layer. Since human dataset is too large, which may take long time, here, we use fly dataset as an example: 
``` 
# load necessary datasets if not loaded: 

data(x_train_fly) 
data(y_train_fly) 
data(x_test_fly) 
data(y_test_fly) 

# producing position frequency matrix for each kernel: 

DeepCirCode_getMotifs(x_train_fly, y_train_fly, x_test_fly, y_test_fly) 

# This step may take tens of minutes since training epochs = 80
# When finished, the "DeepCirCode_Position_Frequency_Matrix.txt" file will be produced in **your current working directory** 
``` 









