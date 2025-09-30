#ifndef NEWTON_H
#define NEWTON_H

#include <Eigen/Dense>
#include <vector>

// Perform Newton update for the z-subproblem inside ADMM
//
// Parameters:
//   z_init    : starting vector (warm start, z^t)
//   x         : knot locations (x1,...,xn)
//.  n         : number of data points
//   A, B      : integration bounds
//   dk_mat, dk_mat_sq : difference operator matrix D^(k)
//   dk        : difference operator matrix D^(k)
//   alpha     : auxiliary variable α^(t)
//   u         : dual variable u^(t)
//   rho       : ADMM penalty parameter
//   max_iters : maximum number of Newton iterations (default 10)
//   tol       : stopping tolerance for gradient norm (default 1e-4)
//
// Returns:
//   z_new     : updated z after Newton iterations
//
Eigen::VectorXd newton_update(
    const Eigen::VectorXd &z_init,
    const std::vector<double> &x,
    const std::vector<double> &W,
    int n,
    double A,
    double B,
    const Eigen::MatrixXd &dk_mat,
    const Eigen::SparseMatrix<double>& dk_mat_sq,
    const Eigen::VectorXd &alpha,
    const Eigen::VectorXd &u,
    double rho,
    int max_iters = 10,
    double tol = 1e-4
);

#endif // NEWTON_H
