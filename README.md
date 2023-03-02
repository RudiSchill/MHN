# MHN
Basic implementation of [Mutual Hazard Networks](https://academic.oup.com/bioinformatics/article/36/1/241/5524604).

## Usage
```
source("UtilityFunctions.R")
source("ModelConstruction.R")
source("Likelihood.R")
source("RegularizedOptimization.R")

Data <- readRDS(file="data/BreastCancer.rds") 

#Learn an MHN from data with sparsity hyperparameter lambda
Theta <- MHN(Data, lambda=0.1)

#Generate a synthetic dataset from an MHN
SyntheticData <- SampleMHN(Theta, size=1000)
```

## Issues
Mac users may get an error when compiling the C code and should install Xcode Command Line Tools via `xcode-select --install`.
