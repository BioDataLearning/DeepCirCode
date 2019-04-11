# DeepCirCode 
In this study, by applying the state-of-art machine learning techniques, we have developed the first deep learning model, DeepCirCode, to predict back-splicing for human circular RNA (circRNA) formation. DeepCirCode utilizes a convolutional neural network (CNN) with nucleotide sequence as the input. Each kernel in the first CNN layer of DeepCirCode can be regarded as a motif scanner. Relevant features for circRNA formation can be automatically learnt by DeepCirCode. In order to make DeepCirCode-learnt features biologically interpretable, we have implemented a visualization method to represent these features as sequence motifs, some of which match human known motifs involved in RNA splicing, transcription or translation. DeepCirCode has been applied to three species for circRNA back-splicing prediction, including human, mouse and fruit fly. Comparision of DeepCirCode-learnt motifs from these species has also been performed.  

This documentation is part of the supplementary information release for DeepCirCode. For details of this work, users can refer to our paper "**Deep Learning of the Back-splicing Code for Circular RNA Formation**" (J. Wang  and L. Wang, 2019). 

# Requirements 
DeepCirCode is an R package with one necessary function written in Python. To use DeepCirCode, R version >= 3.4 is required. Versions for other packages are suggested. 

- R >= 3.4 (64-bit)

- Python 2.7 (64-bit)

- Keras 2.2.4 in R and Python

- ROCR 1.0-7 in R

- reticulate 1.11.1 in R

- numpy in Python

**NOTE:** before using DeepCirCode, users should verify that the Keras package can be used in R by:
```
library("keras")
model <- keras_model_sequential() # Testing keras, no error for this step!!!!
```

# Installing DeepCirCode 
``` 
# Install and load DeepCirCode: 

devtools::install_github("BioDataLearning/DeepCirCode") 
library("DeepCirCode") 
``` 

# Getting the raw training and test sets 
The raw training and test sets of the three species (human, mouse and fruit fly) have been packaged into DeepCirCode, named as "HumanTrain", "HumanTest", "MouseTrain", "MouseTest", "FlyTrain" and "FlyTest", respectively. All datasets are stored as dataframe of dimension [instance_number,3]. The three columns are "**label**" (0 for negative, 1 for positive),"**RNA_input_seq**" (the raw 200 nucleotide RNA sequences, e.g. AUCC...GACG) and "**encoded_seq**" (one-hot encoded RNA sequence of 800 characters in length, e.g. 100000010001...00010010), which can be checked using: 

``` 
# Load datasets in lazy loading mode:

data(HumanTrain)  # datasets are lazy loaded as variable “HumanTrain” until being used
data(HumanTest) 

# Check datasets:

dim(HumanTrain) 
colnames(HumanTrain) 

HumanTrain$label[1] 
HumanTrain$RNA_input_seq[1] 
HumanTrain$encoded_seq[1] 
``` 

# Converting one-hot encoded sequence into 3D matrix as the input for DeepCirCode 
For convenience, we already packaged the 3D matrix datasets into DeepCirCode. Users can jump to the next section for dataset loading and model testing. This section provides users one option for preparing their own datasets.  

DeepCirCode requires the input dimension to be [instance_number,200,4] for sequences of 200-nt in 4 channels (A,T/U,G,C). **x_train** (features of training set) and **x_test** (features of test set) should be of dimension [instance_number,200,4]; **y_train** (label of training set) and **y_test** (label of test set) should be of dimension [instance_number,2] with **y[ ,1]** representing the negative class and **y[ ,2]** representing the positive class. 

``` 
# Prepare x_train and x_test as 3D matrix:

x_train <- convStringToMatrix(HumanTrain$encoded_seq) 
x_test <- convStringToMatrix(HumanTest$encoded_seq) 

# This process may take several minutes
# x_train[1,2,] will print out the one-hot encoded vector of the 2nd nucleotide of HumanTrain$RNA_input_seq[1] 
``` 

``` 
# Prepare y_train and y_test, use Keras function: 

library(keras) 
y_train <- to_categorical(HumanTrain$label,2) 
y_test <- to_categorical(HumanTest$label,2) 
``` 

# Testing DeepCirCode with test sets 
**Step 1**, loading required packages: 

``` 
library(keras) 
library(ROCR) 
``` 
**Step 2**, loading test sets in 3D matrix: 
  
```  
data(x_test_human)
data(y_test_human) 
``` 
**Step 3**, applying DeepCirCode on the test sets: 

The DeepCirCode models for all three species have been packaged into DeepCirCode named "**DeepCirCode_bestmodel_human.hdf5**", "**DeepCirCode_bestmodel_mouse.hdf5**" and "**DeepCirCode_bestmodel_fly.hdf5**". To perform prediction on the test sets, use function **DeepCirCode_test()** which takes arguments the x_test, y_test, and species(="human" or "mouse" or "fly"): 

``` 
DeepCirCode_test(x_test_human, y_test_human, "human") # don't forget the "" for species argument 

# The DeepCirCode model architecture, confusion matrix, loss, accuracy, ROC curve and ROC_AUC value will be printed. 
``` 

# Motif visualization 
We applied the strategy as described by [DanQ](https://github.com/uci-cbcl/DanQ/issues/9) for motif visualization written in python script, named "**getMotifs.py**". "getMotifs.py" defines a function **DeepCirCode_getMotifs()** which takes input x_train, y_train, x_test, and y_test, and outputs a "**DeepCirCode_Position_Probability_Matrix.txt**" file containing the position frequency matrix for each kernel, which can be directly uploaded to [TOMTOM](http://meme-suite.org/tools/tomtom) for motif comparision and visualization. 

Here, we applied R package "reticulate" to run this python script packaged into DeepCirCode. 

**Step 1**, loading package and "getMotifs.py": 

``` 
library(reticulate) 
path <- system.file("getMotifs.py", package = "DeepCirCode")
source_python(path) 
``` 
**Step 2**, applying DeepCirCode_getMotifs(): 

Motif visualization needs the weights of the first layer of DeepCirCode which is not stored in a saved model. So, DeepCirCode_getMotifs() will train DeepCirCode and then visualize the kernels in the first CNN layer:

``` 
# Load datasets if not loaded: 

data(x_train_human) 
data(y_train_human) 
data(x_test_human) 
data(y_test_human) 
```

```
# Get position probability matrix for each kernel: 

DeepCirCode_getMotifs(x_train_human, y_train_human, x_test_human, y_test_human) 

# This step may take tens of minutes since training epochs = 80
# When finished, the "DeepCirCode_Position_Probability_Matrix.txt" file will be saved in **your current working directory** 
``` 









