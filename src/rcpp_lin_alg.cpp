#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map; using Eigen::MatrixXd; using Eigen::Lower; using Eigen::Upper;

#define DBG(MSG,X) Rprintf("%20s SEXP=<%p>. List=%p\n", MSG, (SEXP)X, &X ) ;

// [[Rcpp::export]]
SEXP crossprodCpp(SEXP A, SEXP B = R_NilValue){
  if (TYPEOF(A) != REALSXP | (!Rf_isNull(B) & TYPEOF(B) != REALSXP)) {
    Rcpp::stop("Non-numeric matrix detected. Please convert integer/strings to numeric");
  }
    const Map<MatrixXd> AA(Rcpp::as<Map<MatrixXd> >(A));
    MatrixXd C;
    //If not null do crossprod. Otherwise do crossprod with self
    if(!Rf_isNull(B)){
      const Map<Eigen::MatrixXd> BB(Rcpp::as<Map<MatrixXd> >(B));
      if(AA.rows() != BB.rows()){
        Rcpp::stop("non-conformable arguments");
      }
      C = AA.adjoint() * BB;
    } else{
      const int n(AA.cols());
      C = MatrixXd(n,n).setZero().selfadjointView<Lower>().rankUpdate(AA.adjoint());
    }
    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP matMult( SEXP A, SEXP B){
  if (TYPEOF(A) != REALSXP | TYPEOF(B) != REALSXP) {
    Rcpp::stop("Non-numeric matrix detected. Please convert integer/strings to numeric");
  }
    const Map<MatrixXd> AA(Rcpp::as<Map<MatrixXd> >(A));
    const Map<MatrixXd> BB(Rcpp::as<Map<MatrixXd> >(B));
    if(AA.cols() != BB.rows()){
      Rcpp::stop("non-conformable arguments");
    }
    
    const Eigen::MatrixXd C(AA * BB);
    return Rcpp::wrap(C);
}


//Solve using the householder QR is a solid all around approach
//https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
// [[Rcpp::export]]
SEXP solveCpp( SEXP A, SEXP B = R_NilValue){
  if (TYPEOF(A) != REALSXP | (!Rf_isNull(B) & TYPEOF(B) != REALSXP)) {
    Rcpp::stop("Non-numeric matrix detected. Please convert integer/strings to numeric");
  }
    const Map<MatrixXd> AA(Rcpp::as<Map<MatrixXd> >(A));
    if(AA.cols() != AA.rows()){
      Rcpp::stop("'A' must be square");
    }
    MatrixXd C;
    if (!Rf_isNull(B)){
      const Map<Eigen::MatrixXd> BB(Rcpp::as<Map<MatrixXd> >(B));
      if(AA.cols() != BB.rows()){
        Rcpp::stop("B' must be compatible with 'A'");
      }
      C = AA.householderQr().solve(BB);
    } else{
      C = AA.inverse();
    }
    return Rcpp::wrap(C);
}


// [[Rcpp::export]]
SEXP cholCpp( SEXP A){
  if (TYPEOF(A) != REALSXP) {
    Rcpp::stop("Non-numeric matrix detected. Please convert integer/strings to numeric");
  }
    Map<MatrixXd> AA(Rcpp::as<Map<MatrixXd> >(A));
    if(AA.cols() != AA.rows()){
      Rcpp::stop("'A' must be square");
    }
    const Eigen::MatrixXd L(AA.llt().matrixU()); //.matrixL() for lower
    return Rcpp::wrap(L);
}

// [[Rcpp::export]]
SEXP pseudoinverseCpp( SEXP A){
  if (TYPEOF(A) != REALSXP) {
    Rcpp::stop("Non-numeric matrix detected. Please convert integer/strings to numeric");
  }
  const Map<MatrixXd> AA(Rcpp::as<Map<MatrixXd> >(A));
  const Eigen::MatrixXd pinv(AA.completeOrthogonalDecomposition().pseudoInverse());
  return Rcpp::wrap(pinv);
}

