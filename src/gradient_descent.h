#ifndef TRENDFILTER_GRADIENT_DESCENT_H
#define TRENDFILTER_GRADIENT_DESCENT_H

#include <Eigen/Core>
#include <Eigen/Sparse>

Eigen::VectorXd gd_update(
    const Eigen::VectorXd &z_init,
    const Eigen::VectorXd &W,
    int n,
    const Eigen::SparseMatrix<double>& dk_mat,
    const Eigen::VectorXd &alpha,
    const Eigen::VectorXd &u,
    double rho,
    int max_iters,
    double tol
);

#endif // TRENDFILTER_GRADIENT_DESCENT_H



