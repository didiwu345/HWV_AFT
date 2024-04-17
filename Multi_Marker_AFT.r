library(MASS)
library(CompQuadForm)
library(survival)
library(matrixStats)
library(matrixcalc)
library(RcppArmadillo)
library(Rcpp)
library(Matrix)
library(lbaft)
library(coxKM)

sourceCpp('utils/aft_comprisk_smalln.cpp')
ibskm = function(X){ return(1 - as.matrix(dist(X, method = "manhattan"))/(2*ncol(X))) }  ### IBS kenrel: Input genotype matrix
cpkm = function(X){ return(tcrossprod(X)) }                                              ### Cross-product kernel: Input genotype matrix
simkm = function(X){ return(exp(-as.matrix(dist(X,method='euclidean')^2)/ncol(X))) }     ### Gaussian kernel: Input covariate matrix
idtkm = function(X){ return(0 + outer(as.numeric(X), as.numeric(X), "==")) }             ### For 1-column vector
plkm = function(X){ return(ncol(X) + tcrossprod(X)) }                                    ### Product linear kernel

sample_data = readRDS('sample_data/sample_data.rds')

G = sample_data[['genotype']]
X = sample_data[['adj_cov']]
status = sample_data[['status']]
survt = sample_data[['survt']]
trunct = sample_data[['trunct']]

source('Multi_Marker_AFT.r')

pval = Multi_Marker_AFT(G=G, H=NULL, GxH=NULL, het_cov=NULL, adj_cov=X, kernel_G='ibs', kernel_H=NULL, kernel_het=NULL, smalln_ind='smalln_adj', trunct=trunct, survt=survt, status=status, BB=500)
pval
