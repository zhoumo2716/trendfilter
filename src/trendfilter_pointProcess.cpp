#include <cmath>
#include <stdexcept>
#include <tuple>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <tvdenoising.h>
#include "utils.h"
#include "linearsystem.h"
#include "kf_utils.h"
#include "newton.h"

// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::COLAMDOrdering<int> Ord;

using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::ArrayXd;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Map;

void admm(int n, int k,
                        const Eigen::VectorXd& y, const NumericVector& xd,
                        const Eigen::ArrayXd& weights,
                        double A, double B,
                        Eigen::Ref<Eigen::VectorXd> theta,
                        Eigen::Ref<Eigen::VectorXd> alpha,
                        Eigen::Ref<Eigen::VectorXd> u,
                        int& iter,
                        double& obj_val,
                        const Eigen::SparseMatrix<double>& dk_mat,
                        const Eigen::SparseMatrix<double>& dk_mat_sq,
                        const Eigen::MatrixXd& denseD,
                        const Eigen::VectorXd& s_seq,
                        double lam,
                        double rho,
                        int max_iter,
                        double tol = 1e-5,
                        int newton_max_iters = 10,                 // default max Newton steps
                        double newton_tol = 1e-4) {


  // Perform ADMM updates
  iter = 0;
  VectorXd theta_old = theta;
  double r_norm = 0.0, s_norm = 0.0;
  std::vector<double> x_std(xd.begin(), xd.end());

    //for (iter = 1; iter < max_iter; iter++) {
    // 1. Theta-update: solve through Newton
    //theta = newton_update(
    //  theta,
    //  xd,
    //  n,
    //  A, B,
    //  dk_mat,
    //  dk_mat_sq,
    //  alpha,
    //  u,
    //  rho,
    //  newton_max_iters,
    //  newton_tol);


    // 2. Alpha-update: solve through TV-denoising
    // alpha_k = argmin lambda *||D1*alpha||_1 + (rho/2) * ||alpha - (Dk * theta - u)||_2^2.
    //alpha = tf_dp(dk_mat * theta - u, lam / rho);
    // 3. U-update: dual update
    //u += (alpha - dk_mat * theta);

    // 4. Check convergence (using a simple norm of the primal residuals)
    //r_norm = (alpha - dk_mat * theta).norm();
    //s_norm = (theta - theta_old).norm();
    //if (r_norm < tol && s_norm < tol) break;
    //theta_old = theta;

    // 5. Iteration
    //iter += 1;

  //}

}


