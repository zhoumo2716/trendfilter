#ifndef PM_MATRIX_H  // Prevents multiple inclusions
#define PM_MATRIX_H

#include <RcppEigen.h>  // Include Rcpp dependency

// Declare the function(s) you implemented in `pm_matrix.cpp`
Eigen::MatrixXd pm_matrix(const Eigen::VectorXd &x, int m1, int m2);

#endif

