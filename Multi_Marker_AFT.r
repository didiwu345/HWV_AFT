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

ibskm = function(X){ return(1 - as.matrix(dist(X, method = "manhattan"))/(2*ncol(X))) }  ### IBS kenrel: Input genotype matrix
cpkm = function(X){ return(tcrossprod(X)) }                                              ### Cross-product kernel: Input genotype matrix
simkm = function(X){ return(exp(-as.matrix(dist(X,method='euclidean')^2)/ncol(X))) }     ### Gaussian kernel: Input covariate matrix
idtkm = function(X){ return(0 + outer(as.numeric(X), as.numeric(X), "==")) }             ### Identity kernel for 1-column vector
plkm = function(X){ return(ncol(X) + tcrossprod(X)) }                                    ### Product linear kernel

### Obtain p-value
Multi_Marker_AFT = function(G, H=NULL, GxH=NULL, het_cov=NULL, adj_cov=NULL, kernel_G, kernel_H=NULL, kernel_het=NULL, smalln_ind='no_smalln_adj',trunct=NULL,survt,status,BB=500){
  # Summary: Function implement HWV_AFT association and interaction test on provided inputs and outputs p-value
  # Input:
  #     G:   Genetic variant matrix
  #     H:   Additional genetic variant matrix (Default: NULL; Only needed in G-G/G-E interaction)
  #     GxH: Interaction terms as input matrix (Default: NULL; Only needed in G-G/G-E interaction)
  #     adj_cov:   Adjustment covariate matrix (Default: NULL)
  #     het_cov:   Covariate matrix indicating genetic heterogeneity
  #     kernel_G: Kernel type   (included: 'ibs', 'cp', 'gaussian', 'identity') for G or Kernel matrix
  #     kernel_H: Kernel type   (included: 'ibs', 'cp', 'gaussian', 'identity') for H or Kernel matrix (Default: NULL; Only needed in G-G/G-E interaction)
  #     kernel_het: Kernel type (included: 'ibs', 'cp', 'gaussian', 'identity') for het_cov or Kernel matrix (Default: NULL; Only needed under genetic heterogeneity)
  #     smalln_ind: 'smalln_adj' or 'no_smalln_adj' (Default)
  #     trunct: truncation time, default is NULL --> no truncation time (OR all truncation time is 0)
  #     survt: event or censor time
  #     status: Delta (0 or 1)
  #     BB: Number of iterations for W, S matrix approximation (Default: 500)
  # Output:
  #     pval: p-value
  
  #-- Genetic similarity
  if(is.matrix(kernel_G)){
    Gsim = kernel_G
  } else if(kernel_G=='ibs'){
    Gsim = ibskm(G)
  } else if(kernel_G=='cp'){
    Gsim = cpkm(G)
  } else {
    return("Warning: kernel_G should be either a matrix, or one of  'ibs','cp'")
  }
  
  #-- Define K matrix
  if(is.null(kernel_het)){
    K = Gsim
  } else if(is.matrix(kernel_het)){
    Fsim = kernel_het
    K = Gsim * (1 + Fsim)
  } else if(kernel_het=='ibs'){
    Fsim = ibskm(het_cov)
    K = Gsim * (1 + Fsim)
  } else if(kernel_het=='cp'){
    Fsim = cpkm(het_cov)
    K = Gsim * (1 + Fsim)
  } else if(kernel_het=='gaussian'){
    Fsim = simkm(het_cov)
    K = Gsim * (1 + Fsim)
  } else if(kernel_het=='identity'){
    Fsim = idtkm(het_cov)
    K = Gsim * (1 + Fsim)
  } else {
    return("Warning: kernel_het should be either a matrix, or one of 'ibs','cp','gaussian','identity'")
  }
  
  c = tryCatch({Krank = rankMatrix(K)[1]}, error=function(e){return(100)}, warning=function(e) {return(100)})
  if(c==100){
    return(c)
  } else {
    Krank = rankMatrix(K)[1]
  }
  
  K_eigen = eigen(K)
  K_evec = as.matrix(K_eigen$vectors)[,1:Krank]
  K_eval = abs(K_eigen$values[1:Krank])
  if(Krank==1){
    E1 = matrix(sqrt(K_eval)*K_evec, nrow=1)
  } else {
    E1 = diag(sqrt(K_eval)) %*% t(K_evec)
  }
  
  if(is.null(adj_cov)){
    beta_est = as.matrix(rep(0,ncol(adj_cov)))
  } else {
    if(is.null(trunct)){trunct = rep(0, nrow(G))} else {trunct = trunct}
    fit = suppressWarnings(fit.GP(Surv(as.matrix(trunct), as.matrix(survt), as.matrix(status)) ~ as.matrix(adj_cov) - 1, vMethod='M1'))
    beta_est = as.matrix(fit$beta)
  }
  
  #-- All time points (ordered)
  tall = as.numeric(sort(unlist(log(survt) - as.matrix(adj_cov) %*% beta_est), decreasing = F))
  
  if(smalln_ind=='no_smalln_adj'){
    
    sourceCpp('utils/aft_comprisk.cpp')
    
    W_mats = mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
    Q_all <- Qall_cpp(E1, W_mats, adj_cov, trunct, survt, tall, status, beta_est)
    B_est = matrix(0, nrow=Krank, ncol=ncol(adj_cov))
    for(i in 1:Krank){
      fitlm = lm(Q_all[,i]/sqrt(nrow(G)) ~ 0 + W_mats)
      B_est[i,] = as.numeric(fitlm$coefficients)
    }
    
    S_mats = mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
    U_all <- Uall_cpp(S_mats, adj_cov, trunct, survt, tall, status, beta_est)
    A_est = matrix(0,nrow=ncol(adj_cov),ncol=ncol(adj_cov))
    for(i in 1:ncol(adj_cov)){
      fitlm = lm(U_all[,i]*sqrt(nrow(G)) ~ 0 + S_mats)
      A_est[i,] = as.numeric(fitlm$coefficients)
    }
    
    #-- Covariance matrix of E1M
    E1M_cov <- R1M_cov_cpp(E1, B_est, A_est, adj_cov, trunct, survt, tall, status, beta_est)
    #-- Eigen decomposition
    Chi_coef = abs(eigen(E1M_cov, symmetric=TRUE, only.values=TRUE)$values)
    #-- Test stat and p-value
    E1M = E1 %*% M_tild_cpp(adj_cov, trunct, survt, tall, status, beta_est)
    test_stat = c(crossprod(E1M))
    pval = davies(test_stat, Chi_coef)$Qq
    
  } else if(smalln_ind=='smalln_adj'){
    
    sourceCpp('utils/aft_comprisk_smalln.cpp')
    
    cause1event = status
    M_tild = M_tild_cpp(adj_cov,trunct,survt,tall,cause1event,beta_est)
    q_star = c((t(M_tild) %*% K %*% M_tild) / crossprod(M_tild))
    Knew = K - q_star*diag(nrow(G))
    Knew_rank = rankMatrix(Knew)[1]
    Knew_eigen = eigen(Knew)
    R1 = as.matrix(Knew_eigen$vectors[,1:Knew_rank])
    Knew_eval = Knew_eigen$values[1:Knew_rank]
    
    W_mats <- mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
    Q_all <- Qall_cpp(t(R1),W_mats,adj_cov,trunct,survt,tall,cause1event,beta_est)
    B_est = matrix(0,nrow=Knew_rank,ncol=ncol(adj_cov))
    for(i in 1:Knew_rank){
      fitlm = lm(Q_all[,i]/sqrt(nrow(G)) ~ 0+W_mats)
      B_m = as.numeric(fitlm$coefficients)
      B_est[i,] = B_m
    }
    
    S_mats = mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
    U_all = Uall_cpp(S_mats,adj_cov,trunct,survt,tall,cause1event,beta_est)
    A_est = matrix(0,nrow=ncol(adj_cov),ncol=ncol(adj_cov))
    for(i in 1:ncol(adj_cov)){
      fitlm = lm(U_all[,i]*sqrt(nrow(G)) ~ 0+S_mats)
      A_p = as.numeric(fitlm$coefficients)
      A_est[i,] = A_p
    }
    
    #-- Covariance matrix of R1M
    R1M_cov = R1M_cov_cpp(t(R1),B_est,A_est,adj_cov,trunct,survt,tall,cause1event,beta_est)
    
    #-- R1M covariance Eigen decomposition
    R1M_cov_eig = eigen(R1M_cov,symmetric=TRUE)
    R1M_cov_evec = as.matrix(R1M_cov_eig$vectors)
    R1M_cov_eval = abs(R1M_cov_eig$values)
    
    #-- Get Chi-square distribution coefficients
    Amat = diag(sqrt(R1M_cov_eval)) %*% t(R1M_cov_evec) %*% diag(Knew_eval) %*% R1M_cov_evec %*% diag(sqrt(R1M_cov_eval))
    Chi_coef = eigen(Amat,symmetric=T,only.values=T)$values
    
    #-- Test stat and p-value
    pval = davies(0, Chi_coef)$Qq
    
  } else {
    return("Warning: smalln_ind should be either 'smalln_adj' or 'no_smalln_adj'")
  }
  return(pval)
}

