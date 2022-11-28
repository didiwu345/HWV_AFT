#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat N_bt_cpp(arma::mat X, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::vec epst = log(survt) - (X * beta);

    arma::mat output(n,tnum,fill::zeros);

    for (int i = 0; i < n; i++){
        int vali = cause1event[i];
        for (int tid = 0; tid < tnum; tid++){
            if (epst[i] <= tall[tid]){
                output(i,tid) = vali;
            }
        }
    }
    return output;
}


// [[Rcpp::export]]
arma::mat dN_bt_cpp(arma::mat X, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::mat N_bt = N_bt_cpp(X,survt,tall,cause1event,beta);

    arma::mat output(n,tnum,fill::zeros);

    output.col(0) = N_bt.col(0);
    for (int tid = 1; tid < tnum; tid++){
        output.col(tid) = N_bt.col(tid) - N_bt.col(tid-1);
    }
    return output;
}


// [[Rcpp::export]]
arma::mat Y_bt_cpp(arma::mat X, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::vec epst = log(survt) - (X * beta);

    arma::mat output(n,tnum,fill::zeros);

    for (int i = 0; i < n; i++){
        for (int tid = 0; tid < tnum; tid++){
            if(epst[i] >= tall[tid]){
                output(i,tid) = 1;
            }
        }
    }

    return output;
}


// [[Rcpp::export]]
arma::mat V1_bt_cpp(arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::vec epsa = log(trunct) - (X * beta);

    arma::mat output(n,tnum,fill::zeros);

    for (int i = 0; i < n; i++){
        for (int tid = 0; tid < tnum; tid++){
            if(epsa[i] <= tall[tid]){
                output(i,tid) = 1;
            }
        }
    }
    return output;
}


// [[Rcpp::export]]
arma::vec colSum(arma::mat XX){
    int p = XX.n_cols;
    arma::vec output(p,fill::zeros);
    for(int j = 0; j < p; j++){
        output[j] = sum(XX.col(j));
    }
    return output;
}


// [[Rcpp::export]]
arma::vec M_tild_cpp(arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;

    arma::mat dN_bt = dN_bt_cpp(X, survt, tall, cause1event, beta);
    arma::mat V1_bt = V1_bt_cpp(X, trunct, survt, tall, cause1event, beta);
    arma::mat Y_bt = Y_bt_cpp(X, survt, tall, cause1event, beta);
    arma::mat VY_bt = V1_bt % Y_bt;
    arma::vec sum_dN_sum_VY = colSum(dN_bt)/colSum(VY_bt);

    arma::vec output(n,fill::zeros);

    for (int i=0; i<n; i++){
        output[i] = sum(dN_bt.row(i) - VY_bt.row(i)%sum_dN_sum_VY.as_row());
    }

    return output;

}


// [[Rcpp::export]]
arma::mat Qall_cpp(arma::mat R1, arma::mat W_mats, arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;
    int BB = W_mats.n_rows;
    int m = R1.n_rows;

    arma::vec Q = R1 * M_tild_cpp(X,trunct,survt,tall,cause1event,beta);

    arma::mat output(BB,m,fill::zeros);

    for (int i = 0; i < BB; i++){
        arma::vec newbeta = beta + W_mats.row(i).as_col() / sqrt(n);
        arma::vec newmat = R1 * M_tild_cpp(X,trunct,survt,tall,cause1event,newbeta) - Q;
        output.row(i) = newmat.as_row();
    }
    return output;
}



// [[Rcpp::export]]
arma::vec Ubeta_cpp(arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;
    int p = X.n_cols;

    arma::mat dN_bt = dN_bt_cpp(X, survt, tall, cause1event, beta);
    arma::mat V1_bt = V1_bt_cpp(X, trunct, survt, tall, cause1event, beta);
    arma::mat Y_bt = Y_bt_cpp(X, survt, tall, cause1event, beta);
    arma::mat VY_bt = V1_bt % Y_bt;
    arma::vec sum_VY = colSum(VY_bt);
    arma::mat sum_XVY = X.t() * VY_bt;

    arma::vec output(p,fill::zeros);

    for (int i = 0; i < n; i++){
        arma::vec Xi = X.row(i).as_col();
        for (int tid = 0; tid < tnum; tid++){
            if (sum_VY[tid] > 0){
                output += (Xi - sum_XVY.col(tid)/sum_VY[tid]) * dN_bt(i,tid);
            }
        }
    }
    return output / n;
}



// [[Rcpp::export]]
arma::mat Uall_cpp(arma::mat S_mats, arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;
    int BB = S_mats.n_rows;
    int p = X.n_cols;

    arma::mat output(BB,p,fill::zeros);

    for (int i = 0; i < BB; i++){
        arma::vec newbeta = beta + S_mats.row(i).as_col() / sqrt(n);
        arma::vec newmat = Ubeta_cpp(X, trunct, survt, tall, cause1event, newbeta);
        output.row(i) = newmat.as_row();
    }
    return output;
}



// [[Rcpp::export]]
arma::mat Vmat_cpp(arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::mat V1_bt = V1_bt_cpp(X, trunct, survt, tall, cause1event, beta);
    arma::mat Y_bt = Y_bt_cpp(X, survt, tall, cause1event, beta);
    arma::mat VY_bt = V1_bt % Y_bt;
    arma::vec sum_VY = colSum(VY_bt);

    arma::mat output(n,tnum,fill::zeros);

    for (int i = 0; i < n; i++){
        for (int tid = 0; tid < tnum; tid++){
            if (sum_VY[tid] > 0){
                output(i,tid) = VY_bt(i,tid)/sum_VY[tid];
            }
        }
    }
    return output;    
}



// [[Rcpp::export]]
arma::mat dLambda_cpp(arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::mat dN_bt = dN_bt_cpp(X, survt, tall, cause1event, beta);
    arma::mat V1_bt = V1_bt_cpp(X, trunct, survt, tall, cause1event, beta);
    arma::mat Y_bt = Y_bt_cpp(X, survt, tall, cause1event, beta);
    arma::mat VY_bt = V1_bt % Y_bt;
    arma::vec sum_VY = colSum(VY_bt);
    arma::vec sum_dN = colSum(dN_bt);

    arma::mat output(n,tnum,fill::zeros);

    for (int i = 0; i < n; i++){
        for (int tid = 0; tid < tnum; tid++){
            if (sum_VY[tid]>0){
                output(i,tid) = VY_bt(i,tid) * sum_dN[tid] / sum_VY[tid];
            }
        }
    }
    return output;
}



// [[Rcpp::export]]
arma::mat R1M_cov_cpp(arma::mat R1, arma::mat B_est, arma::mat A_est, arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec cause1event, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;
    int p = A_est.n_rows;

    arma::mat dN_bt = dN_bt_cpp(X, survt, tall, cause1event, beta);
    arma::mat Vmat = Vmat_cpp(X, trunct, survt, tall, cause1event, beta);
    arma::mat dLambda = dLambda_cpp(X, trunct, survt, tall, cause1event, beta);
    arma::vec sum_dN = colSum(dN_bt);

    arma::mat output2(n,n,fill::zeros);

    for (int tid = 0; tid < tnum; tid++){
        output2 += diagmat(dLambda.col(tid)) - (Vmat.col(tid) * Vmat.col(tid).as_row()) * sum_dN[tid];
    }

    arma::mat output1 = R1 - B_est * solve(A_est,eye(p,p)) * trans(X);

    return output1 * output2 * trans(output1);
}


