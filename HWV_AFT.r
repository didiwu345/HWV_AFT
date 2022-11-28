
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

seed_list = read.table("~/Documents/Genetics_aft/random_seeds.txt")[,1]

source("~/Documents/Genetics_aft/SNP.R")
real.genome = readRDS('~/Documents/Genetics_aft/dat.rds')

sourceCpp('~/Documents/Genetics_aft/aft_comprisk.cpp')


### Data simulation function
data_simulation = function(size, p, valtype, adj, trunc){
    """
    # size: sample size
    # p: SNP size
    # valtype: 't1e' or 'power'
    # adj: w/ or w/o adjustment ('adj' or 'noadj')
    # trunc: w/ or w/o left truncation ('trunc' or 'notrunc')
    """
    
    I = diag(size)
    ones = as.matrix(rep(1,size))
    
    #-- SNP-set (Generate SNPs from 1000-Genom)
    real.geno = snpdata(n=size, p=p, real.genome)
    G = real.geno$G
    maf = real.geno$mafs
    
    #-- Covariates
    if(adj=='adj'){
        X1 = rbinom(size,1,prob=0.5)
        X2 = runif(size,0,2)
    } else {
        X1 = rep(0,size)
        X2 = rep(0,size)
    }
    X = as.matrix(cbind(X1,X2))
    eps = rnorm(size)
    
    beta1 = as.matrix(rep(0.035, p))
    alpha1 = as.matrix(rep(0.015, ncol(X)))
    beta2 = as.matrix(rep(0.1, p))
    alpha2 = as.matrix(rep(0.03, ncol(X)))
    
    #-- gene/covariate effect
    beta1G = c(G %*% beta1)
    alpha1X = c(X %*% alpha1)
    beta2G = c(G %*% beta2)
    alpha2X = c(X %*% alpha2)
    
    # cause1 sum(Gbeta + Xalpha)
    if(valtype == 't1e'){   sum_xbeta1 =          alpha1X }
    if(valtype == 'power'){ sum_xbeta1 = beta1G + alpha1X }
    # cause2 sum(Gbeta + Xalpha)
    sum_xbeta2 = beta2G + alpha2X
    
    # cause-specific hazard
    csh1 = function(valtype,t,id){ return((t[id]^2*exp(-3*sum_xbeta1[id])+t[id]*exp(-2*sum_xbeta1[id]))) }
    csh2 = function(t,id){    return((t[id]^2*exp(-3*sum_xbeta2[id])+t[id]*exp(-2*sum_xbeta2[id]))) }
    
    #-- Generate Sruvival times
    func = function(t){
        expsum1 = (exp(-3*sum_xbeta1[id]) + exp(-3*sum_xbeta2[id]))
        expsum2 = (exp(-2*sum_xbeta1[id]) + exp(-2*sum_xbeta2[id]))
        return(1/3 * expsum1 * t^3 + 1/2 * expsum2 * t^2 + log(cdf))
    }
    
    survt = NULL
    trunct = NULL
    id = 1
    while(id<=size){
        if(trunc=='trunc'){
          newtrunct = runif(1,0,1)
        } else {
          newtrunct = 0
        }
        cdf = runif(1,0,1)
        newsurvt = uniroot(func, c(0,100000))$root
        while(newsurvt %in% survt){
          cdf  = runif(1,0,1)
          newsurvt = uniroot(func, c(0,100000))$root
        }
        if(newsurvt > newtrunct){
          trunct = c(trunct,newtrunct)
          survt = c(survt,newsurvt)
          id = id + 1
        }
    }
    
    #-- determine interest OR competing event type
    event_prob = apply(as.matrix(c(1:size)), 1, function(x){ csh1(valtype,survt,x)/(csh1(valtype,survt,x)+csh2(survt,x)) })
    cause1_idx = NULL
    for(k in 1:length(survt)){ cause1_idx = c(cause1_idx, rbinom(1,1,event_prob[k])) }
    
    #-- times for event of interest
    times_event1 = survt[cause1_idx==1]
    #-- generate censor times
    if(trunc=='trunc'){
        censor = trunct[cause1_idx==1] + rexp(sum(cause1_idx), 0.1)
    } else {
        censor = rexp(sum(cause1_idx), 0.4)
    }
    #-- surv_time = min(times, censor)
    survt[cause1_idx==1] = apply(cbind(times_event1,censor), 1, min)
    #-- right-censored event of interest
    status = cause1_idx
    status[status==1] = (times_event1 < censor) + 0
    censor_rate = sum(status==0)/length(status)
    
    output = list()
    output[['G']] = G
    output[['X']] = X
    output[['censor_rate']] = sum(status)/length(status)
    output[['status']] = status
    output[['survt']] = survt
    output[['trunct']] = trunct
    
    return(output)
       
}


ibskm = function(X){
    return(1 - as.matrix(dist(X, method = "manhattan"))/(2*ncol(X))) 
}  ### IBS kenrel: Input genotype matrix

cpkm = function(X){
    return(tcrossprod(X))
}  ### Cross-product kernel: Input genotype matrix

### Obtain p-value
HWV_AFT = function(G, H=NULL, GxH=NULL, het_cov=NULL, adj_cov=NULL, kernel_G, kernel_H=NULL, kernel_het=NULL, smalln_ind='no_smalln_adj', BB=400){
    """
    Summary: Function implement HWV_AFT association and interaction test on provided inputs and outputs p-value
    Input:
        G:   Genetic variant matrix
        H:   Additional genetic variant matrix (Default: NULL; Only needed in G-G/G-E interaction)
        GxH: Interaction terms as input matrix (Default: NULL; Only needed in G-G/G-E interaction)
        adj_cov:   Adjustment covariate matrix (Default: NULL)
        het_cov:   Covariate matrix indicating genetic heterogeneity
        kernel_G: Kernel type   ('ibs', 'cp') for G or Kernel matrix
        kernel_H: Kernel type   ('ibs', 'cp') for H or Kernel matrix (Default: NULL; Only needed in G-G/G-E interaction)
        kernel_het: Kernel type ('ibs', 'cp') for het_cov or Kernel matrix (Default: NULL; Only needed under genetic heterogeneity)
        smalln_ind: 'smalln_adj' or 'no_smalln_adj' (Default)
        BB: Number of iterations for W, S matrix approximation (Default: 400)
    Output:
        pval: p-value
    """
    
    if(smalln_ind == 'smalln_adj')
    
    #-- Genetic similarity
    if(is.matrix(kernel_G)){
        Gsim = kernel_G
    } else if(kernel_G=='ibs'){
        Gsim = ibskm(G)
    } else if(kernel_G=='cp'){
        Gsim = cpkm(G)
    } else {
        return("Warning: kernel_G should be either a matrix, 'ibs' or 'cp'")
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
    } else {
        return("Warning: kernel_het should be either a matrix, 'ibs' or 'cp'")
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
        beta_est = as.matrix(rep(0,ncol(X)))
    } else if(is.matrix(adj_cov)){
        fit = suppressWarnings(fit.GP(Surv(trunct, survt, status) ~ adj_cov - 1, vMethod='M1'))
        beta_est = as.matrix(fit$beta)
    }
    
    #-- All time points (ordered)
    tall = sort(c(log(survt) - adj_cov %*% beta_est), decreasing = F)
    
    if(smalln_ind='no_smalln_adj'){
        
        W_mats = mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
        Q_all <- Qall_cpp(E1, W_mats, adj_cov, trunct, survt, tall, status, beta_est)  #time:0.5s
        B_est = matrix(0, nrow=Krank, ncol=ncol(adj_cov))
        for(i in 1:Krank){
            fitlm = lm(Q_all[,i]/sqrt(nrow(G)) ~ 0 + W_mats)
            B_est[i,] = as.numeric(fitlm$coefficients)
        }
        
        S_mats = mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
        U_all <- Uall_cpp(S_mats, adj_cov, trunct, survt, tall, status, beta_est)   #time:0.5s
        A_est = matrix(0,nrow=ncol(adj_cov),ncol=ncol(adj_cov))
        for(i in 1:ncol(X)){
            fitlm = lm(U_all[,i]*sqrt(nrow(G)) ~ 0 + S_mats)
            A_est[i,] = as.numeric(fitlm$coefficients)
        }
        
        #-- Covariance matrix of E1M
        E1M_cov <- R1M_cov_cpp(E1, B_est, A_est, adj_cov, trunct, survt, tall, status, beta_est)   #time:0.2s
        #-- Eigen decomposition
        Chi_coef = abs(eigen(E1M_cov, symmetric=TRUE, only.values=TRUE)$values)
        #-- Test stat and p-value
        E1M = E1 %*% M_tild_cpp(adj_cov, trunct, survt, tall, status, beta_est)
        test_stat = c(crossprod(E1M))
        pval = davies(test_stat, Chi_coef)$Qq
        
    } else if(smalln_ind='smalln_adj'){
        
        W_mats = mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
        Q_all = Qall_cpp(t(R1),W_mats,adj_cov,trunct,survt,tall,cause1event,beta_est)  # 17s
        B_est = matrix(0,nrow=Knew_rank,ncol=ncol(adj_cov))
        for(i in 1:Knew_rank){
            fitlm = lm(Q_all[,i]/sqrt(nrow(G)) ~ 0+W_mats)
            B_m = as.numeric(fitlm$coefficients)
            B_est[i,] = B_m
        }
        
        S_mats = mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
        U_all = Uall_cpp(S_mats,adj_cov,trunct,survt,tall,cause1event,beta_est)  # 5.5s
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



