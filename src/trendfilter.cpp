#include <cmath>
#include <stdexcept>
#include <tuple>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "utils.h"
extern "C"{
  #include "tf_dp.h"
}

// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::COLAMDOrdering<int> Ord;

using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::ArrayXd;
using Eigen::VectorXd;
using Eigen::MatrixXd;

/**
 * Implementation of linear system solve using Sparse QR, tridiagonal, or
 * other updates (e.g., Kalman Filter-based).
 *
 * Encapsulating all approaches to solving the linear system in the theta
 * update allows a single ADMM routine to be written (rather than
 * several separate routines, one for each solver).
 *
 * Additional solvers to include in the future: Kalman-Filter based...
 */
class LinearSystem {
  SparseQR<SparseMatrix<double>, Ord> qr;
  VectorXd a, b, c, cp;
  public:
    void compute(const SparseMatrix<double>&, bool);
    std::tuple<VectorXd,int> solve(const VectorXd&, bool);
};
void LinearSystem::compute(const SparseMatrix<double>& A, bool tridiag) {
  if (tridiag) {
    std::tie(a, b, c) = extract_tridiag(A);
    cp = tridiag_forward(a, b, c);
  } else {
    // Setting the pivot threshold to be negative forces SparseQR to be
    // maximally conservative in dropping columns, which is important when
    // the Gram matrix is ill-conditioned (which often is the case for
    // unequally spaced inputs.
    qr.setPivotThreshold(-1.0);
    qr.compute(A);
  }
}
std::tuple<VectorXd,int> LinearSystem::solve(const VectorXd& v, bool tridiag) {
  if (tridiag) {
    return std::make_tuple(tridiag_backsolve(a, b, cp, v), 0);
  } else {
    VectorXd sol = qr.solve(v);
    int info = int(qr.info());
    return std::make_tuple(sol, info);
  }
}

double tf_gauss_loss(const VectorXd& y,
                     const VectorXd& theta,
                     const ArrayXd& weights) {
  return 0.5 * ((y - theta).array() * weights.sqrt()).matrix().squaredNorm();
}

double tf_penalty(const VectorXd& theta, const NumericVector& xd, double lam, int k) {
  Eigen::VectorXd Dv = Dkv(theta, k + 1, xd, true);
  return lam * Dv.lpNorm<1>();
}


// Gaussian only
double tf_objective(const VectorXd& y, const VectorXd& theta,
                    const NumericVector& xd,
                    const ArrayXd& weights, double lam, int k) {
  return tf_gauss_loss(y, theta, weights) + tf_penalty(theta, xd, lam, k);
}

// Matches initialization from package glmgen
// Not used in these ADMM updates (projection via Legendre polynomials is
// more stable, but left here as a reference implementation).
std::tuple<VectorXd,int> init_theta_nullspace(const VectorXd& y, const
    SparseMatrix<double>& penalty_mat) {
  SparseQR<SparseMatrix<double>, Ord> qr;
  qr.compute(penalty_mat.transpose());
  // Project y onto rowspace of C^{k+1}
  VectorXd beta = qr.solve(y);
  // Subtract residual off rowspace-projection to obtain nullspace projection.
  return std::make_tuple(y - (penalty_mat.transpose()) * beta, int(qr.info()));
}

VectorXd init_u(const VectorXd& residual, const NumericVector& xd, int k,
    const ArrayXd& weights) {
  // (Dk)^T u = y - \theta, from stationarity
  SparseQR<SparseMatrix<double>, Ord> qr;
  SparseMatrix<double> dk_mat = get_dk_mat(k, xd, false);
  qr.compute(dk_mat.transpose());
  return qr.solve((residual.array()*weights).matrix());
}

void admm_single_lambda(
    int n,
    const Eigen::VectorXd& y,
    const NumericVector& xd,
    const Eigen::ArrayXd& weights,
    int k,
    Eigen::VectorXd theta,
    Eigen::VectorXd alpha,
    Eigen::VectorXd u,
    int& iter,
    double& obj_val,
    const Eigen::SparseMatrix<double>& dk_mat_sq,
    double lam,
    int max_iter,
    double rho,
    double tol = 1e-5,
    bool tridiag = false) {
  // Initialize internals
  VectorXd tmp(n-k);
  VectorXd Dth_tmp(alpha.size());
  VectorXd alpha_old(alpha);
  VectorXd wy = (y.array()*weights).matrix();
  double rr, ss;

  // Form Gram matrix and set up linear system for theta update
  SparseMatrix<double> A = rho * dk_mat_sq;
  A.diagonal().array() += weights;
  A.makeCompressed();

  LinearSystem linear_system;
  // Technically, can form one SparseQR object, analyze the pattern once,
  // and then re-use it.
  // So call analyzePattern once, and then factorize repeatedly.
  // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseQR.html#aba8ae81fd3d4ce9139eccb6b7a0256b2
  linear_system.compute(A, tridiag);

  // Perform ADMM updates
  int computation_info;
  iter = 0;
  double best_objective = tf_objective(y, theta, xd, weights, lam, k);
  VectorXd best_theta = theta;
  for (iter = 1; iter < max_iter; iter++) {
    if (iter % 1000 == 0) Rcpp::checkUserInterrupt(); // check if killed

    // theta update

    std::tie(theta, computation_info) = linear_system.solve(
        wy + rho * Dktv(alpha + u, k, xd), tridiag);
    // if (computation_info > 1) {
    //  std::cerr << "Eigen Sparse QR solve returned nonzero exit status.\n";
    // }
    Dth_tmp = Dkv(theta, k, xd);
    tmp = Dth_tmp - u;
    // alpha update
    tf_dp(n-k, tmp.data(), lam/rho, alpha.data());
    // u update
    u += alpha - Dth_tmp;
    double cur_objective = tf_objective(y, theta, xd, weights, lam, k);
    if (cur_objective < best_objective) {
      best_objective = cur_objective;
      best_theta = theta;
    }

    // Check for convergence
    rr = (Dth_tmp - alpha).norm() / alpha.size();
    ss = rho * Dktv(alpha - alpha_old, k, xd).norm() / dk_mat_sq.rows();
    alpha_old = alpha;
    if (rr < tol && ss < tol) break;
  }
  // For checking the criteria. Generally, the dual (ss) is much slower.
  // Rcpp::Rcout << iter << ": rr = " << rr << " ss = " << ss << std::endl;
}

/*
 * Tech debt: because tf_dp is implemented in C and expects y.data() (and
 * weights.data(), cannot qualify y as a const even though we never intend to
 * modify it.
 */
// [[Rcpp::export]]
Rcpp::List admm_lambda_seq(
    NumericVector x,
    Eigen::VectorXd y,
    Eigen::ArrayXd weights,
    int k,
    Eigen::VectorXd lambda,
    int nlambda = 50,
    double lambda_max = -1.0,
    double lambda_min = -1.0,
    double lambda_min_ratio = 1e-5,
    int max_iter = 200,
    double rho_scale = 1.0,
    double tol = 1e-5,
    bool tridiag = false) {

  int n = x.size();

  if (lambda[0] < tol / 100 && lambda_max <= 0) {
    lambda_max = get_lambda_max(x, y, weights, k);
  }
  get_lambda_seq(lambda, lambda_max, lambda_min, lambda_min_ratio, nlambda);


  Eigen::MatrixXd theta(n, nlambda);
  Rcpp::NumericVector objective_val(nlambda);
  Rcpp::IntegerVector iters(nlambda);

  // Use DP solution for k=0.
  if (k == 0) {
    for (int i = 0; i < nlambda; i++) {
      tf_dp_weight(n, y.data(), weights.data(), lambda[i],
                   theta.col(i).data());
      objective_val[i] = tf_objective(y, theta.col(i), x, weights, lambda[i], k);
    }
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("theta") = theta,
      Rcpp::Named("lambda") = lambda,
      Rcpp::Named("tf_objective") = objective_val,
      Rcpp::Named("iters") = iters
    );
    return out;
  }


  // Initialize difference matrices and other helper objects
  SparseMatrix<double> dk_mat = get_dk_mat(k, x, false);
  SparseMatrix<double> dk_mat_sq = dk_mat.transpose() * dk_mat;
  Eigen::MatrixXd alpha(n-k, nlambda);

  // Initialize ADMM variables
  // Project onto Legendre polynomials to initialize for largest lambda.
  theta.col(0) = project_polynomials(x, y, weights, k);
  alpha.col(0) = Dkv(theta.col(0), k, x);
  VectorXd u = init_u((theta.col(0) - y)/(lambda[0]*rho_scale), x, k, weights);


  for (int i = 0; i < nlambda; i++) {
    Rcpp::checkUserInterrupt();
    admm_single_lambda(n, y, x, weights, k,
      theta.col(i), alpha.col(i), u, // return vals
      iters[i], objective_val[i],
      dk_mat_sq, lambda[i], max_iter, lambda[i]*rho_scale,
      tol, tridiag);
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("theta") = theta,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("lambda") = lambda,
    Rcpp::Named("tf_objective") = objective_val,
    Rcpp::Named("iters") = iters
  );
  return out;
}


Rcpp::List admm_single_lambda_with_tracking(NumericVector x,
    Eigen::VectorXd& y, const Eigen::ArrayXd& weights, int k,
    double lam, int max_iter, double rho,
    bool tridiag=false) {

  int n = x.size();

  // Initialize difference matrices and other helper objects
  SparseMatrix<double> penalty_mat = get_penalty_mat(k+1, x);
  SparseMatrix<double> dk_mat = get_dk_mat(k, x, false);
  VectorXd wy = (y.array()*weights).matrix();

  // Initialize theta, alpha, u
  VectorXd theta = project_polynomials(x, y, weights, k);
  VectorXd alpha = dk_mat * theta;
  VectorXd u = init_u((theta-y)/rho, x, k, weights);
  VectorXd tmp(n-k);

  // Form Gram matrix and set up linear system for theta update
  SparseMatrix<double> A = rho * (dk_mat.transpose() * dk_mat);
  A.diagonal().array() += weights;
  A.makeCompressed();

  if (tridiag & (k != 1)) {
    throw std::invalid_argument("`tridiag` can only be used with k=1.");
  }
  LinearSystem linear_system;
  linear_system.compute(A, tridiag);

  // Perform ADMM updates
  int computation_info;
  MatrixXd theta_mat(n, max_iter);
  VectorXd objective_vec(max_iter);
  int iter = 0;
  theta_mat.col(iter) = theta;
  objective_vec[iter] = tf_objective(y, theta, x, weights, lam, k);
  VectorXd best_theta = theta;
  double best_objective = objective_vec[iter];
  for (iter = 1; iter < max_iter; iter++) {
    // theta update
    std::tie(theta, computation_info) = linear_system.solve(
        wy + rho * (dk_mat.transpose() * (alpha + u)),
        tridiag);
    // if (computation_info > 1) {
    //   std::cerr << "Eigen Sparse QR solve returned nonzero exit status.\n";
    // }
    tmp = dk_mat*theta - u;
    // alpha update
    tf_dp(n-k, tmp.data(), lam/rho, alpha.data());
    // u update
    u += alpha - dk_mat*theta;
    objective_vec[iter] = tf_objective(y, theta, x, weights, lam, k);
    theta_mat.col(iter) = theta;
    if (objective_vec[iter] < best_objective) {
      best_objective = objective_vec[iter];
      best_theta = theta;
    }
  }
  // TODO: Implement stopping criterion
  Rcpp::List return_list = Rcpp::List::create(
      Rcpp::Named("theta")=best_theta,
      Rcpp::Named("objective")=objective_vec,
      Rcpp::Named("theta_mat")=theta_mat,
      Rcpp::Named("computation_info")=computation_info);
  return return_list;
}
