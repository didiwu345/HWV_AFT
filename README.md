# HWV_AFT


'HWV_AFT.r' is the main function that implements HWV_AFT association and interaction test.

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

'sample_data' folder contains essential sample datasets to conduct an example association test, without considering genetic heterogeneity, with covariate adjustmne and left truncation.

'util' folder contains utility C++ functions 'aft_comprisk.cpp' and 'aft_comprisk_smalln.cpp' that will be used in main function 'HWV_AFT.r'.



