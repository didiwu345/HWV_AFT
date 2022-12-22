# Multi_Marker_AFT

## Introduction
`Multi_Marker_AFT.r` is the main function that implements multi-marker association and interaction test based on the AFT model.

To successfully implement the test, you will need following R package dependencies:
* MASS
* CompQuadForm
* survival
* matrixStats
* matrixcalc
* RcppArmadillo
* Rcpp
* Matrix
* lbaft
* coxKM

`sample_data` folder contains essential sample datasets to conduct an example association test, without considering genetic heterogeneity, with covariate adjustment and left truncation.

`utils` folder contains utility C++ functions `aft_comprisk.cpp` and `aft_comprisk_smalln.cpp` that will be used in main function `Multi_Marker_AFT.r`.

`Missing_genotype_imputation` folder has a sample script to perform missing genotype imputation with IMPUTE2 software. The best practices using IMPUTE2 software can be found below:
[Getting started with IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#getting_started)

## Example of Multi_Marker_AFT method implementation
### Step1: Clone the repo
```bash
cd ~
git clone https://github.com/didiwu345/Multi_Marker_AFT.git
cd Multi_Marker_AFT
```

### Step2: Sample implementation
This code example implements the `Multi_Marker_AFT` method in detecting genetic association effects, with covariate adjustment and left truncation, without considering confounding effect and genetic heterogeneity. 
The sample genotype dataset has sample size `N=500`, and SNP-set size `p=10`. Since `N >> p`, small sample correction is turned off. `IBS` kernel is used to measure genetic similarity.

```bash
suppressMessages(library(MASS))
suppressMessages(library(CompQuadForm))
suppressMessages(library(survival))
suppressMessages(library(matrixStats))
suppressMessages(library(matrixcalc))
suppressMessages(library(RcppArmadillo))
suppressMessages(library(Rcpp))
suppressMessages(library(Matrix))
suppressMessages(library(lbaft))
suppressMessages(library(coxKM))

sourceCpp('utils/aft_comprisk.cpp')
ibskm = function(X){ return(1 - as.matrix(dist(X, method = "manhattan"))/(2*ncol(X))) }  ### IBS kenrel: Input genotype matrix
cpkm = function(X){ return(tcrossprod(X)) }                                              ### Cross-product kernel: Input genotype matrix
simkm = function(X){ return(exp(-as.matrix(dist(X,method='euclidean')^2)/ncol(X))) }     ### Gaussian kernel: Input covariate matrix
idtkm = function(X){ return(0 + outer(as.numeric(X), as.numeric(X), "==")) }             ### For 1-column vector
plkm = function(X){ return(ncol(X) + tcrossprod(X)) }                                    ### Product linear kernel

G = read.csv('sample_data/genotype.csv')
X = read.csv('sample_data/adj_cov.csv')
status = read.csv('sample_data/status.csv')
survt = read.csv('sample_data/survt.csv')
trunct = read.csv('sample_data/trunct.csv')

pval = Multi_Marker_AFT(G=G, H=NULL, GxH=NULL, het_cov=NULL, adj_cov=X, kernel_G='ibs', kernel_H=NULL, kernel_het=NULL, smalln_ind='no_smalln_adj', trunct=trunct, survt=survt, status=status, BB=500)
```

