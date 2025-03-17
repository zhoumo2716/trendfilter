#include <Rcpp.h>
#include <RcppEigen.h>
#include "pm_matrix.h"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;


Eigen::MatrixXd compute_P_matrix(const Eigen::VectorXd &x_target, const Eigen::VectorXd &x_support) {
  int m = x_target.size();
  int n = x_support.size();
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(m, n);

  for (int i = 0; i < m; i++) {
    P(i, 0) = 1.0;
    for (int j = 1; j < n; j++) {
      P(i, j) = P(i, j - 1) * (x_target[i] - x_support[j - 1]);
    }
  }
  return P;
}

Eigen::MatrixXd compute_C_withoutTheta(const Eigen::VectorXd &x_support) {
  int n = x_support.size();
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n, n);

  for (int t = 0; t < n; t++) {
    for (int w = 0; w <= t; w++) {
      double prod_denom = 1.0;
      for (int u = 0; u <= t; u++) {
        if (u != w) {
          prod_denom *= (x_support[w] - x_support[u]);
        }
      }
      C(t, w) = 1.0 / prod_denom;
    }
  }
  return C;
}

Eigen::MatrixXd compute_A_matrix(const Eigen::VectorXd &x_target, const Eigen::VectorXd &x_support) {
  Eigen::MatrixXd P = compute_P_matrix(x_target, x_support);
  Eigen::MatrixXd C_withoutTheta = compute_C_withoutTheta(x_support);

  int m = P.rows();
  int n = C_withoutTheta.cols();

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, n);

  A = P * C_withoutTheta;

  // for (int i = 0; i < m; i++) {
  //   for (int j = 0; j < n; j++) {
  //     A(i, j) = 0;
  //     for (int k = 0; k < n; k++) {
  //       A(i, j) += P(i, k) * C_withoutTheta(k, j);
  //     }
  //   }
  // }
  return A;
}

// [[Rcpp::export]]
Eigen::MatrixXd pm_matrix(const Eigen::VectorXd &x, int m1, int m2) {
  int n = x.size();
  if ((m1 + m2) > n) Rcpp::stop("Error: m1 + m2 must be <= length of x");

  // Define target and support points
  Eigen::VectorXd x_target_left = x.head(m1);
  Eigen::VectorXd x_support_left = x.segment(m1, m1);
  Eigen::VectorXd x_target_right = x.tail(m2);
  Eigen::VectorXd x_support_right = x.segment(n-2*m2, m2);

  // Compute A matrices for left and right
  Eigen::MatrixXd A_left = compute_A_matrix(x_target_left, x_support_left);
  Eigen::MatrixXd A_right = compute_A_matrix(x_target_right, x_support_right);
  Eigen::MatrixXd A_middle = Eigen::MatrixXd::Identity(n - m1 - m2, n - m1 - m2); // Identity matrix for middle


  // Initialize full A matrix
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n - m1 - m2);

  // Fill left, middle, and right portions
  A.block(0, 0, A_left.rows(), A_left.cols()) = A_left;
  A.block(m1, 0, A_middle.rows(), A_middle.cols()) = A_middle;
  A.block(n - m2, n - m1 - 2*m2, A_right.rows(), A_right.cols()) = A_right;

  // for (int i = 0; i < A_left.nrow(); i++) {
  //   for (int j = 0; j < A_left.ncol(); j++) {
  //     A(i, j) = A_left(i, j);
  //   }
  // }
  //
  // for (int i = 0; i < A_middle.nrow(); i++) {
  //   for (int j = 0; j < A_middle.ncol(); j++) {
  //     A(m1 + i, j) = A_middle(i, j);
  //   }
  // }
  //
  // for (int i = 0; i < A_right.nrow(); i++) {
  //   for (int j = 0; j < A_right.ncol(); j++) {
  //     A(n - m2 + i, n - m1 - 2*m2 + j) = A_right(i, j);
  //   }
  // }

  return A;
}
