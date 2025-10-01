#ifndef INTEGRAL_WEIGHTS_H
#define INTEGRAL_WEIGHTS_H

#include <Eigen/Dense>

// [[Rcpp::depends(RcppEigen)]]

Eigen::VectorXd compute_integral_weights(
    const Eigen::VectorXd &x,
    double A,
    double B,
    int n
);

#endif // INTEGRAL_WEIGHTS_H

