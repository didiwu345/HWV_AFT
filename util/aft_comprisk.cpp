#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat N_bt_cpp(arma::mat X, arma::vec survt, arma::vec tall, arma::vec status, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::vec epst = log(survt) - (X * beta);

    arma::mat output(n,tnum,fill::zeros);

    for (int i = 0; i < n; i++){
        int vali = status[i];
        for (int tid = 0; tid < tnum; tid++){
            if (epst[i] <= tall[tid]){
                output(i,tid) = vali;
            }
        }
    }
    return output;
}


// [[Rcpp::export]]
arma::mat dN_bt_cpp(arma::mat X, arma::vec survt, arma::vec tall, arma::vec status, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::mat N_bt = N_bt_cpp(X,survt,tall,status,beta);

    arma::mat output(n,tnum,fill::zeros);

    output.col(0) = N_bt.col(0);
    for (int tid = 1; tid < tnum; tid++){
        output.col(tid) = N_bt.col(tid) - N_bt.col(tid-1);
    }
    return output;
}


// [[Rcpp::export]]
arma::mat V1Y_bt_cpp(arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::vec epst = log(survt) - (X * beta);
    arma::vec epsa = log(trunct) - (X * beta);

    arma::mat output(n,tnum,fill::zeros);

    for (int i = 0; i < n; i++){
        for (int tid = 0; tid < tnum; tid++){
            if(epst[i] >= tall[tid] && epsa[i] <= tall[tid]){
                output(i,tid) = 1;
            }
        }
    }

    return output;
}


// [[Rcpp::export]]
arma::vec M_tild_cpp(arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec status, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::mat dN_bt = dN_bt_cpp(X, survt, tall, status, beta);
    arma::mat VY_bt = V1Y_bt_cpp(X, trunct, survt, tall, beta);
    arma::vec onevec = ones(tnum);

    return dN_bt*onevec - VY_bt*((trans(dN_bt)*onevec)/(trans(VY_bt)*onevec));

}


// [[Rcpp::export]]
arma::mat Qall_cpp(arma::mat R1, arma::mat W_mats, arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec status, arma::vec beta){

    int n = X.n_rows;
    int BB = W_mats.n_rows;
    int m = R1.n_rows;

    arma::vec Q = R1 * M_tild_cpp(X,trunct,survt,tall,status,beta);

    arma::mat output(BB,m,fill::zeros);

    for (int i = 0; i < BB; i++){
        arma::vec newbeta = beta + W_mats.row(i).as_col() / sqrt(n);
        arma::vec newmat = R1 * M_tild_cpp(X,trunct,survt,tall,status,newbeta) - Q;
        output.row(i) = newmat.as_row();
    }
    return output;
}


// [[Rcpp::export]]
arma::vec Ubeta_cpp(arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec status, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;
    int p = X.n_cols;

    arma::mat dN_bt = dN_bt_cpp(X, survt, tall, status, beta);
    arma::mat VY_bt = V1Y_bt_cpp(X, trunct, survt, tall, beta);
    arma::vec onevec = ones(tnum);

    arma::vec output1 = (trans(X)*dN_bt)*onevec;
    arma::mat XVY = trans(X)*VY_bt;
    arma::vec colsum_VY = trans(VY_bt)*onevec;

    arma::vec output2(p,fill::zeros);
    for (int q=0; q<p; q++){
        output2[q] = as_scalar(trans(onevec)*(dN_bt*(XVY.row(q).as_col()/colsum_VY)));
    }
    
    return (output1-output2)/n;
}



// [[Rcpp::export]]
arma::mat Uall_cpp(arma::mat S_mats, arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec status, arma::vec beta){

    int n = X.n_rows;
    int BB = S_mats.n_rows;
    int p = X.n_cols;

    arma::mat output(BB,p,fill::zeros);

    for (int i = 0; i < BB; i++){
        arma::vec newbeta = beta + S_mats.row(i).as_col() / sqrt(n);
        arma::vec newmat = Ubeta_cpp(X, trunct, survt, tall, status, newbeta);
        output.row(i) = newmat.as_row();
    }
    return output;
}



// [[Rcpp::export]]
arma::mat Vmat_cpp(arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec status, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::mat VY_bt = V1Y_bt_cpp(X, trunct, survt, tall, beta);
    arma::vec onevec = ones(tnum);
    arma::vec sum_VY = trans(VY_bt)*onevec;

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
arma::mat dLambda_cpp(arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec status, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;

    arma::mat dN_bt = dN_bt_cpp(X, survt, tall, status, beta);
    arma::mat VY_bt = V1Y_bt_cpp(X, trunct, survt, tall, beta);
    arma::vec onevec = ones(tnum);
    arma::vec sum_VY = trans(VY_bt)*onevec;
    arma::vec sum_dN = trans(dN_bt)*onevec;

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
arma::mat R1M_cov_cpp(arma::mat R1, arma::mat B_est, arma::mat A_est, arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec status, arma::vec beta){

    int n = X.n_rows;
    int tnum = tall.n_elem;
    int p = A_est.n_rows;

    arma::mat dN_bt = dN_bt_cpp(X, survt, tall, status, beta);
    arma::mat Vmat = Vmat_cpp(X, trunct, survt, tall, status, beta);
    arma::mat dLambda = dLambda_cpp(X, trunct, survt, tall, status, beta);
    arma::vec onevec = ones(tnum);
    arma::vec sum_dN = trans(dN_bt)*onevec;

    arma::mat output2(n,n,fill::zeros);

    for (int tid = 0; tid < tnum; tid++){
        output2 += diagmat(dLambda.col(tid)) - (Vmat.col(tid) * Vmat.col(tid).as_row()) * sum_dN[tid];
    }

    arma::mat output1 = R1 - B_est * solve(A_est,eye(p,p)) * trans(X);

    return output1 * output2 * trans(output1);
}


// // [[Rcpp::export]]
// arma::mat R1M_cov_cppnew(arma::mat R1, arma::mat B_est, arma::mat A_est, arma::mat X, arma::vec trunct, arma::vec survt, arma::vec tall, arma::vec status, arma::vec beta){

//     int n = X.n_rows;
//     int tnum = tall.n_elem;
//     int p = A_est.n_rows;

//     arma::mat dN_bt = dN_bt_cpp(X, survt, tall, status, beta);
//     arma::mat Vmat = Vmat_cpp(X, trunct, survt, tall, status, beta);
//     arma::mat dLambda = dLambda_cpp(X, trunct, survt, tall, status, beta);
//     arma::vec onevec = ones(tnum);
//     arma::vec sum_dN = trans(dN_bt)*onevec;
//     arma::vec colsum_dLambda = trans(dLambda)*onevec;

//     arma::mat output2(n,n,fill::zeros);

//     for (int tid=0; tid<tnum; tid++){
//         arma::mat mata = dLambda.col(tid)*trans(Vmat.col(tid));
//         output2 += diagmat(dLambda.col(tid)) - mata - trans(mata) + Vmat.col(tid)*trans(Vmat.col(tid))*colsum_dLambda[tid]; 
//     }

//     arma::mat output1 = R1 - B_est * solve(A_est,eye(p,p)) * trans(X);

//     return output1 * output2 * trans(output1);
// }


