#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "utils.h"
#include "kf_utils.h"

// [[Rcpp::depends(RcppEigen)]]

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::ArrayXd;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Map;

// in-place computation of Ptemp given P and A
Eigen::MatrixXd computePtemp(Eigen::MatrixXd A, Eigen::MatrixXd P) {
  int k = P.rows();
  Eigen::MatrixXd temp = A.row(0) * P;
  Eigen::MatrixXd var = A.row(0) * temp.transpose();

  // in-place block replacement starting from the last component:
  P.block(1, 1, k - 1, k - 1).reverse() = P.block(0, 0, k - 1, k - 1).reverse();
  P(0, 0) = var(0, 0);
  temp.conservativeResize(1, k - 1);  // drop the last component
  P.block(0, 1, 1, k - 1) = temp;
  P.block(1, 0, k - 1, 1) = temp.transpose();
  return P;
}

// [[Rcpp::export]]
Eigen::MatrixXd smat_to_mat(const Eigen::SparseMatrix<double>& sparseMat, int k, bool equal_spaced) {
  int rows = sparseMat.rows(); // n-k
  Eigen::MatrixXd denseMat(rows, k + 1);
  std::vector<double> rowNonzeros;
  // Iterate over nonzero coefficients in each row of the sparse matrix
  for (int i = 0; i < sparseMat.outerSize(); ++i) {
    std::vector<double> rowNonzeros;
    for (SparseMatrix<double>::InnerIterator it(sparseMat, i); it; ++it) 
      rowNonzeros.push_back(it.value());
    int m = rowNonzeros.size();
    for (int j = 0; j < m; j++) {
      if (i < rows) 
        denseMat(i - j, j) = rowNonzeros[m - 1 - j];
      else  
        denseMat(rows - 1 - j, j + i - rows + 1) = rowNonzeros[m - 1 - j];
    }
    if (equal_spaced && i == k) {
      denseMat.conservativeResize(1, k + 1);
      return denseMat;
    } 
  }
  return denseMat;
}

void configure_denseD(Rcpp::NumericVector x, Eigen::MatrixXd& denseD, Eigen::VectorXd& s_seq, Eigen::SparseMatrix<double>& dk_mat, int k, bool equal_space) {
    int n = x.size();
    // construct dense D mat with only nonzero entries from sparse D mat:
    denseD = smat_to_mat(dk_mat, k, equal_space);

    // ideally, construct dense denseD directly:
    // denseD = b_mat(k, x, VectorXi::LinSpaced(n - k, 0, n - k - 1));
    // re-construct denseD in which each row saves values for T.row(0) per iterate: 
    if (equal_space) {  // if x is equally spaced, the dense D matrix only contains the nonzero band.
      s_seq = denseD.block(0, k, 1, 1);
      denseD.conservativeResize(1, k);
    } else {
      //s_seq = VectorXd::Ones(n - k);
      s_seq = denseD.block(0, k, n - k, 1);
      denseD.conservativeResize(n - k, k);
    }
    // reverse order
    MatrixXd firstRow(1, k);
    int m = denseD.rows();
    for (int i = 0; i < m; i++) {
      firstRow = -denseD.row(i) / s_seq(i);
      std::reverse(firstRow.data(), firstRow.data() + k);
      denseD.row(i) = firstRow;
    }
}

// [[Rcpp::export]]
Rcpp::List configure_denseD_test(Rcpp::NumericVector x, int k) {
    SparseMatrix<double> dk_mat = get_dk_mat(k, x, false);
    int n = x.size();
    bool equal_space = is_equal_space(x, 0.1);
    MatrixXd denseD = equal_space ? MatrixXd::Zero(1, k) : MatrixXd::Zero(n, k);
    VectorXd s_seq = equal_space ? VectorXd::Zero(1) : VectorXd::Zero(n);
  
    configure_denseD(x, denseD, s_seq, dk_mat, k, equal_space);
    
    return Rcpp::List::create(Rcpp::Named("s_seq") = s_seq, Rcpp::Named("D_mat") = denseD);
}

void f1step(double y, double c, double Z, double H, const Eigen::MatrixXd& A, 
  double RQR, Eigen::VectorXd& a, Eigen::MatrixXd& P, double& vt, double& Ft,
  Eigen::VectorXd& Kt) {
  VectorXd a_temp = A * a;
  a_temp(0) += c;
  MatrixXd Ptemp = computePtemp(A, P);  // A * P * A.transpose();
  Ptemp(0, 0) += RQR;

  vt = y - Z * a_temp(0);
  Ft = pow(Z, 2) * Ptemp(0, 0) + H;
  Kt = Ptemp.col(0) * Z;
  a = a_temp + Kt * vt / Ft;
  P = Ptemp - Kt * Z * Ptemp.row(0) / Ft;

  // symmetrize
  Ptemp = P;
  Ptemp += P.transpose();
  P = Ptemp / 2;

  // Some entries of P can be _really_ small in magnitude, set them to 0
  // but maintain symmetric / posdef
  for (int i = 0; i < P.rows(); i++) {
    for (int j = 0; j <= i; j++) {
      if (abs(P(i, j)) < 1e-30) {
        if (i == j) {
          P.row(i).setZero();
          P.col(i).setZero();
        } else {
          P(i, j) = 0;
          P(j, i) = 0;
        }
      }
    }
  }
}

void df1step(double y, double Z, double H, const Eigen::MatrixXd& A, double RQR,
  Eigen::VectorXd& a, Eigen::MatrixXd& P, Eigen::MatrixXd& Pinf, int& rankp,
  double& vt, double& Ft, double& Finf, Eigen::VectorXd& Kt, Eigen::VectorXd& Kinf) {
  double tol = Eigen::NumTraits<double>::epsilon();
  tol = std::sqrt(tol);
  int k = a.size();
  MatrixXd Ptemp(k, k);

  a = A * a;
  P = computePtemp(A, P);  // A * P * A.transpose();
  P(0, 0) += RQR;
  Pinf = computePtemp(A, Pinf);  // A * Pinf * A.transpose();

  vt = y - Z * a(0);
  Kt = P.col(0) * Z;
  Ft = Z * Kt(0) + H;
  Kinf = Pinf.col(0) * Z;
  Finf = Z * Kinf(0);

  if (Finf > tol) {  // should always happen
    a += vt * Kinf / Finf;
    P += Ft * Kinf * Kinf.transpose() / pow(Finf, 2);
    P -= Kt * Kinf.transpose() / Finf + Kinf * Kt.transpose() / Finf;
    Pinf -= Kinf * Kinf.transpose() / Finf;
    rankp--;
  } else {  // should never happen
    Finf = 0;
    if (Ft > tol) {
      a += vt * Kt / Ft;
      P -= Kt * Kt.transpose() / Ft;
    }
  }
  if (Ft < tol)
    Ft = 0;

  // symmetrize
  Ptemp = P;
  Ptemp += P.transpose();
  P = Ptemp / 2;
  Ptemp = Pinf;
  Ptemp += Pinf.transpose();
  Pinf = Ptemp / 2;

  // Fix possible negative definiteness, should never happen
  for (int i = 0; i < k; i++) {
    if (P(i, i) < 0) {
      P.row(i).setZero();
      P.col(i).setZero();
    }
    if (Pinf(i, i) < 0) {
      Pinf.row(i).setZero();
      Pinf.col(i).setZero();
    }
  }
}
