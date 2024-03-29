# metNet

## metNet Overview

This repository is the code base and model depo for the classification of cancers of unknown primary using methylation data. This 
model extends from 450k data to the classification of the InfiniumEpicArray (850K)

## Installation
We utilized multiple programming languages (Python, and R) to construct learning models 
and to perform biological analyses. We have provided a snakemake pipeline to ingest and 
analyze the data. The features will be filtered for maximally informative subset of features
differentiating  tumor class.

**Docker Installation**
#### Not active yet TODO:
The docker image for this project can be pulled from the online Docker Hub [incomplete build] 
or can be built using the Dockerfile included in the base directory of this project.
To pull the image from the Docker Hub repo, run the following command:
```
docker pull mskaro1/metNet
```

To build the image using the Dockerfile, run the following command in the base directory of this project:
```
docker build --tag mskaro1/metNet .
```

**Recommneded: Docker Approach**
```
docker run --rm -it -v <output-directory> mskaro1/metNet
```
Note: `<output-directory>` should be replaced with the path of a directory on the user's local machine, and it is where the outputs of the demo will be stored.
Note: This command should be run in the base directory of the project, and `<output-directory>` should be replaced with the path to a directory for the outputs to be stored.
  
**Outputs**
```
├── <output-directory>/
│   ├── model/
│   ├── oversampled-datasets/
|   ├── feature-selected-datasets/
|   ├── classification-results/
```
- model: Saved h5. 
- oversampled-datasets: The training and testing data generated from the binary datasets. The training data uses synthetic data generated by the SMOTE algorithm, while the testing data uses only real TCGA data.
- important-features: The top ±5000 features (i.e. probes) of each training dataset ranked by their RF regression score.
- feature-selected-datasets: The training and testing datasets that only contain the top 1000 selected features.
- classification-results: Directory contains the classification results of our Random Forest model on the feature-selected datasets.

## General Usage:
---------
* metNet, when run without the :code:

```
usage: metNet.py -i INPUT -o OUTPUT -m MODEL -c classifications
ML pipeline to classify tissue organ of origin for cancers of unknown primary using illumina 450k and minfinum epic array data

required arguments:
          -i INPUT, --input
Path to input directory containing the methylation data
          -o OUTPUT, --output
Path to output directory to save results
          -m MODEL, --model
Path to CUP pretrained model
          -c CLASS_DIR, --classifications
Path to classifications visualization data.
          -s STATE --state
raw or normalized Binary variable indicating if the data in the input directory is raw idat 
data or normalized beta values. 
```                              
    
The `-h` flag to understand all available options. 

Additionally, each component of the pipeline can be called individually from the command line. 

**Note:** For those seeking to use the docker image to interact with our framework, run the following command to gain access to the shell of the docker image:
```
docker run --rm -it --entrypoint="" mskaro1/metNet bash
```
## Reviewing:

Docker image in prep ! 
TODO: [] Check out our wiki for implementation actions!
Thanks!
