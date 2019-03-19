# SAVERX

R package for transfer learning of scRNA-seq denoising. Take a look at our free [SAVER-X web-server](https://singlecell.wharton.upenn.edu/saver-x/) for the transfer learning online computation! We also encourage you to read our [pre-print manucript](https://www.biorxiv.org/content/10.1101/457879v2) for more information. You can also refer to our earlier denoising method [SAVER](http://github.com/mohuangx/SAVER).

## Installation

First, install the supporting Python package sctransfer. See the source code of the package [here](http://github.com/jingshuw/sctransfer)

```
pip install sctransfer
```

Next, open R and install the R package SAVERX
```
library(devtools)
install_github("jingshuw/SAVERX")
```

## Basic Usage

Our current pre-trained models can be downloaded [here](https://www.dropbox.com/sh/4u22cfuswcfcwvu/AAC6CgsO7dvQSNInTF0wWMDva?dl=0)

Our input can be '.txt', '.csv' or '.rds' file. The '.rds' file can store either a matrix or a sparse matrix of class 'dgCMatrix'. 
As a toy example, you may download one of our demo datasets on the web server, shekhar_downsampled.csv, the down-sampled mouse retina data from [here](https://www.dropbox.com/sh/kctbw41kdh6jmnb/AAAO5Icu97Ep6uoWFdHRKIcMa?dl=0). Say you have saved the file in a folder './testdata/'. As SAVER-X will generate intemediate files in the same folder as the input dataset during the computation, please make sure that only one SAVER-X task is running on the folder. 

### SAVER-X without pretraining

```
library(SAVERX)
saverx("./testdata/shekhar_downsampled.csv")
```

### SAVER-X with a pretrained model
For the demo dataset, we have a pre-trained model for the mouse retina, please download the file, mouse_Retina.hdf5, and you may save it in './mouse_retina.hdf5'
```
library(SAVERX)
saverx("./testdata/shekhar_downsampled.csv", data.species = "Mouse", 
use.pretrain = T, pretrained.weights.file = "./mouse_retina.hdf5", model.species = "Mouse")
```

For both cases, you will find the final results in './testdata/shekhar_downsampled_denoised.rds'. When dealing with large datasets, you can set 'is.large.data = T' to reduce RAM.

