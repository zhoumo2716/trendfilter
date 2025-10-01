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
#include "integral_weights.h"
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



void admm(int n,
          int k,
          const Eigen::VectorXd& xd,
          const Eigen::VectorXd& W,
          double A, double B,
          Eigen::VectorXd& theta,
          Eigen::VectorXd& alpha,
          Eigen::VectorXd& u,
          const Eigen::SparseMatrix<double>& dk_mat,
          const Eigen::SparseMatrix<double>& dk_mat_sq,
          double lam,
          double rho,
          int& iter,
          int max_iter = 100,
          double tol = 1e-5,
          int newton_max_iters = 50,  // default max Newton steps
          double newton_tol = 1e-4) {


  // Perform ADMM updates
  Eigen::VectorXd theta_old = theta;
  double r_norm = 0.0, s_norm = 0.0;

  for (iter = 1; iter < max_iter; iter++) {
    // 1. Theta-update: solve through Newton
    theta = newton_update(
      theta,
      W,
      n,
      A, B,
      dk_mat,
      dk_mat_sq,
      alpha,
      u,
      rho,
      newton_max_iters,
      newton_tol);


    // 2. Alpha-update: solve through TV-denoising
    // alpha_k = argmin lambda *||D1*alpha||_1 + (rho/2) * ||alpha - (Dk * theta - u)||_2^2.
    alpha = tf_dp(dk_mat * theta - u, lam / rho);
    // 3. U-update: dual update
    u += (alpha - dk_mat * theta);

    // 4. Check convergence (using a simple norm of the primal residuals)
    r_norm = (alpha - dk_mat * theta).norm();
    s_norm = (theta - theta_old).norm();
    if (r_norm < tol && s_norm < tol) break;
    theta_old = theta;

  }

}


// [[Rcpp::export]]
Rcpp::List trendfilter_pointProcess(NumericVector x,
                     int k,
                     double A, double B,
                     double lambda = 1,
                     double rho_scale = 1,
                     int max_iter = 100,
                     double tol = 1e-5,
                     int newton_max_iters = 50,
                     double newton_tol = 1e-5) {

  // Project onto Legendre polynomials to initialize for largest lambda.
  // Convert R -> Eigen
  VectorXd xd = Rcpp::as<VectorXd>(x);
  int n = static_cast<int>(xd.size());
  int dim = n+2;

  Rcpp::NumericVector x_aug(dim);
  x_aug[0] = A;
  std::copy(xd.begin(), xd.end(), x_aug.begin() + 1);
  x_aug[dim - 1] = B;

  // Initialize difference matrices and other helper objects
  Eigen::SparseMatrix<double> dk_mat = get_dk_mat(k, x_aug, false);
  Eigen::SparseMatrix<double> dk_mat_sq = dk_mat.transpose() * dk_mat;

  Eigen::VectorXd y = Eigen::VectorXd::Ones(dim);
  Eigen::VectorXd weights = Eigen::VectorXd::Ones(dim);
  Eigen::VectorXd theta = project_polynomials(x_aug, y, weights, k);
  Eigen::VectorXd alpha = Dkv(theta, k, x);

  // For the dual, a common simple choice is to initialize it to zero:
  Eigen::VectorXd u = Eigen::VectorXd::Zero(alpha.size());


  // Integral weights
  Eigen::VectorXd W = compute_integral_weights(xd, A, B, n);


  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////

  // Rcpp::Rcout << "[debug] sizes:"
  //             << "  n=" << n
  //             << "  dim=" << n+2
  //             << "  theta=" << theta.size()
  //             << "  alpha=" << alpha.size()
  //             << "  u="     << u.size()
  //             << "  W="     << W.size()
  //             << "\n";
  //
  // Rcpp::Rcout << "[debug] dk_mat dims:    "
  //             << dk_mat.rows() << " x " << dk_mat.cols() << "\n";
  // Rcpp::Rcout << "[debug] dk_mat_sq dims: "
  //             << dk_mat_sq.rows() << " x " << dk_mat_sq.cols() << "\n";
  //
  // R_FlushConsole();
  // R_ProcessEvents();
  //Rcpp::stop("Debug abort before ADMM: see printed dimensions above.");

  ///////////////////////////////////////////////////////////////////////////////////

  int iter = 0;
  admm(n, k, xd, W, A, B, theta, alpha, u,
       dk_mat, dk_mat_sq, lambda, lambda*rho_scale,
       iter, max_iter, tol, newton_max_iters, newton_tol);

  return Rcpp::List::create(
    Rcpp::Named("theta") = theta,
    Rcpp::Named("iter") = iter
  );
}


