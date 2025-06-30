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
#include "pm_matrix.h"
#include "softThreshold.h"

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

void admm_single_lambda_D1(int n, int k,
                        const Eigen::VectorXd& y, const NumericVector& xd,
                        Eigen::Ref<Eigen::VectorXd> theta,
                        Eigen::Ref<Eigen::VectorXd> alpha1, Eigen::Ref<Eigen::VectorXd> alphak,
                        Eigen::Ref<Eigen::VectorXd> u1, Eigen::Ref<Eigen::VectorXd> uk,
                        const Eigen::SparseMatrix<double>& dk_mat,
                        const Eigen::SparseMatrix<double>& d1_mat,
                        double rho,
                        double lambda1, double lambdak,
                        int max_iter,
                        double tol,
                        int &iters) {

  // Form the LHS matrix in theta update: A = I + rho*(Dk^T * Dk) + rho*(D1^T * D1)
  MatrixXd lhs = MatrixXd::Identity(n, n) + rho * (dk_mat.transpose() * dk_mat) + rho * (d1_mat.transpose() * d1_mat);


  VectorXd theta_old = theta;
  double r_norm = 0.0, s_norm = 0.0;

  Eigen::LDLT<MatrixXd> solver;
  solver.compute(lhs);

  // ADMM iterations
  for (int iter = 0; iter < max_iter; iter++) {
    // 1. Theta-update: solve for theta (our beta variable)
    //    (I + rho*(Dk^T Dk) + rho*(D1^T D1)) * theta = y + rho*Dk^T(alpha_k + u_k) + rho*D1^T(alpha_1 + u_1)
    VectorXd rhs = y + rho * (dk_mat.transpose() * (alphak + uk)) + rho * (d1_mat.transpose() * (alpha1 + u1));
    theta = solver.solve(rhs);

    // 2. Alpha_k update (for the D^(k+1) penalty)
    //    Solve: alpha_k = argmin lambda_k1 * ||D1*alpha_k||_1 + (rho/2) * ||alpha_k - (Dk * theta - u_k)||_2^2.
    //    1D fussed Lasso
    VectorXd zk = dk_mat * theta - uk;
    alphak = tf_dp(zk, lambdak / rho);

    // 3. Alpha_1 update (for the D^(1) penalty)
    //    Solve: alpha_1 = argmin lambda_1 * ||alpha_1||_1 + (rho/2) * ||alpha_1 - (D1 * theta - u_1)||_2^2.
    //    SoftThresholding solver
    VectorXd z1 = d1_mat * theta - u1;
    alpha1 = softThreshold(z1, lambda1 / rho);

    // 4. Dual updates
    uk += (alphak - dk_mat * theta);
    u1 += (alpha1 - d1_mat * theta);

    // 5. Check convergence (using a simple norm of the primal residuals)
    r_norm = std::max((alphak - dk_mat * theta).norm(), (alpha1 - d1_mat * theta).norm());
    s_norm = (theta - theta_old).norm();
    if (r_norm < tol && s_norm < tol) break;
    theta_old = theta;

    // 6. Iteration
    iters += 1;
  }
}

// [[Rcpp::export]]
Rcpp::List trendfilter_D1(int k, double lambda1_scalar, double lambdak, const Eigen::VectorXd& y, const NumericVector& x, double rho_scale = 1.0) {
  SparseMatrix<double> dk_mat = get_dk_mat(k, x, false);
  SparseMatrix<double> d1_mat = get_dk_mat(1, x, false);

  int n = y.size();
  Eigen::ArrayXd weights = Eigen::ArrayXd::Ones(n);

  // Project onto Legendre polynomials to initialize for largest lambda.

  Eigen::VectorXd theta(n);

  theta = project_polynomials(x, y, weights, k);

  Eigen::VectorXd alphak = Dkv(theta, k, x);
  Eigen::VectorXd alpha1 = Dkv(theta, 1, x);

  // For the dual, a common simple choice is to initialize it to zero:
  Eigen::VectorXd uk = Eigen::VectorXd::Zero(alphak.size());
  Eigen::VectorXd u1 = Eigen::VectorXd::Zero(alpha1.size());


  double lambda1 = lambda1_scalar * lambdak;
  double rho = rho_scale*(lambda1+lambdak)/2;
  double tol = 1e-5;
  int max_iter = 1000;
  int iters = 0;

  admm_single_lambda_D1(n, k,
                        y, x,
                        theta,
                        alpha1, alphak,
                        u1, uk,
                        dk_mat,
                        d1_mat,
                        rho,
                        lambda1, lambdak,
                        max_iter,
                        tol,
                        iters);

// Rcpp::List out = Rcpp::List::create(
//     Rcpp::Named("theta") = theta,
//     Rcpp::Named("alpha1") = alpha1,
//     Rcpp::Named("alphak") = alphak,
//     Rcpp::Named("iters") = iters);
//

Rcpp::NumericMatrix theta_mat(theta.size(), 1);
std::copy(theta.begin(), theta.end(), theta_mat.begin());

Rcpp::NumericVector lambda_vec(1);
lambda_vec[0] = lambdak;

  List out = List::create(
    _["x"]      = x,
    _["theta"]  = theta_mat,
    _["lambda"] = lambda_vec,
    _["lambda1_scalar"] = lambda1_scalar,
    _["iters"] = iters
  );



  out.attr("class") = CharacterVector::create("trendfilter_D1", "trendfilter");
  return out;

}

