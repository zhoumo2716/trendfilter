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
using Eigen::Map;

/**
 * Implementation of linear system solve using Sparse QR, tridiagonal, or
 * other updates (e.g., Kalman Filter-based).
 *
 * Encapsulating all approaches to solving the linear system in the theta
 * update allows a single ADMM routine to be written (rather than
 * several separate routines, one for each solver).
 *
 * @param solver = 0 (`trigiag`), 1 (`sparse_qr`), 2 (`kf`)
 */
class LinearSystem {
  public:
    void construct(const VectorXd&, const Eigen::ArrayXd&, int, double, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&, int);
    void compute(int);
    std::tuple<VectorXd,int> solve(const VectorXd&, int, NumericVector, double, int);
    void kf_init(const VectorXd&, const ArrayXd&, int , double, const SparseMatrix<double>&);
    void kf_iter(const Eigen::VectorXd&);
  private:
    // tridiag
    VectorXd a, b, c, cp;
    VectorXd wy;
    // sparse_qr
    SparseQR<SparseMatrix<double>, Ord> qr;
    SparseMatrix<double> A;
    // kf: 
    int d, rankp;
    double vt_b, Ft_b, Finf_b;
    VectorXd wsqrt, yw, RQR, a1, vt, Ft, Finf, Kt_b, Kinf_b, s_seq, r,  r1, rtmp, sol;
    MatrixXd T, at, P1, Pt, P1inf, Pinf, Kt, Kinf, L0, L1, Ptemp, Dseq;
};

void LinearSystem::construct(const Eigen::VectorXd& y, const Eigen::ArrayXd& weights, int k, double rho, const Eigen::SparseMatrix<double>& dk_mat, const Eigen::SparseMatrix<double>& dk_mat_sq, int solver) {
  switch(solver) {
    case 0: {
      wy = (y.array() * weights).matrix();
      // Form Gram matrix and set up linear system for theta update
      A = rho * dk_mat_sq;
      A.diagonal().array() += weights;
      A.makeCompressed();
      break;
    }
    case 1: {
      wy = (y.array() * weights).matrix();
      // Form Gram matrix and set up linear system for theta update
      A = rho * dk_mat_sq;
      A.diagonal().array() += weights;
      A.makeCompressed();
      break;
    }
    case 2: {
      LinearSystem::kf_init(y, weights, k, rho, dk_mat);
      break;
    }
  } 
}

void LinearSystem::compute(int solver) {
  switch(solver) {
    case 0: {
      std::tie(a, b, c) = extract_tridiag(A);
      cp = tridiag_forward(a, b, c);
      break;
    }
    case 1: {
      // Setting the pivot threshold to be negative forces SparseQR to be
      // maximally conservative in dropping columns, which is important when
      // the Gram matrix is ill-conditioned (which often is the case for
      // unequally spaced inputs.
      qr.setPivotThreshold(-1.0);
      qr.compute(A);
      break;
    }
    case 2:
      break;
  }
}
std::tuple<VectorXd,int> LinearSystem::solve(const Eigen::VectorXd& adj_mean, int k, Rcpp::NumericVector x, double rho, int solver) {
  int info = 0;
  switch(solver) {
    case 0: {
      VectorXd v = wy + rho * Dktv(adj_mean, k, x); // c = alpha + u
      sol = tridiag_backsolve(a, b, cp, v);
      break;
    }
    case 1: {
      VectorXd v = wy + rho * Dktv(adj_mean, k, x);
      sol = qr.solve(v);
      info = int(qr.info());
      break;
    }
    case 2: {
      LinearSystem::kf_iter(adj_mean);
      break;
    }
  }
  return std::make_tuple(sol, info);
}

// [[Rcpp::export]]
Eigen::VectorXd linear_single_solve_test(int linear_solver, const Eigen::VectorXd y, const Eigen::ArrayXd weights, const Rcpp::NumericVector x, double rho, const Eigen::VectorXd adj_mean) {
  int n = y.size();
  int k = n - adj_mean.size();
  SparseMatrix<double> dk_mat = get_dk_mat(k, x, false);
  SparseMatrix<double> dk_mat_sq = dk_mat.transpose() * dk_mat;
  
  VectorXd sol = VectorXd::Zero(n);
  LinearSystem linear_system;
  int info = 0;
  linear_system.construct(y, weights, k, rho, dk_mat, dk_mat_sq, linear_solver);
  linear_system.compute(linear_solver);
  std::tie(sol, info) = linear_system.solve(adj_mean, k, x, rho, linear_solver);
  return sol;
}

LinearSystem linear_system;

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
    Eigen::Ref<Eigen::VectorXd> theta,
    Eigen::Ref<Eigen::VectorXd> alpha,
    Eigen::Ref<Eigen::VectorXd> u,
    int& iter,
    double& obj_val,
    const Eigen::SparseMatrix<double>& dk_mat,
    const Eigen::SparseMatrix<double>& dk_mat_sq,
    double lam,
    int max_iter,
    double rho,
    double tol = 1e-5,
    int linear_solver = 2) {
  // Initialize internals
  VectorXd tmp(n-k);
  VectorXd Dth_tmp(alpha.size());
  VectorXd alpha_old(alpha);
  VectorXd wy = (y.array()*weights).matrix();
  double rr, ss;

  // LinearSystem linear_system;
  // Technically, can form one SparseQR object, analyze the pattern once,
  // and then re-use it.
  // So call analyzePattern once, and then factorize repeatedly.
  // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseQR.html#aba8ae81fd3d4ce9139eccb6b7a0256b2
  linear_system.construct(y, weights, k, rho, dk_mat, dk_mat_sq, linear_solver);
  linear_system.compute(linear_solver);

  // Perform ADMM updates
  int computation_info;
  iter = 0;
  // double best_objective = tf_objective(y, theta, xd, weights, lam, k);
  // VectorXd best_theta = theta;
  for (iter = 1; iter < max_iter; iter++) {
    if (iter % 1000 == 0) Rcpp::checkUserInterrupt(); // check if killed

    // theta update
    std::tie(theta, computation_info) = linear_system.solve(
        alpha + u, k, xd, rho, linear_solver);
    // if (computation_info > 1) {
    //  std::cerr << "Eigen Sparse QR solve returned nonzero exit status.\n";
    // }
    Dth_tmp = Dkv(theta, k, xd);
    tmp = Dth_tmp - u;
    // alpha update
    tf_dp(n-k, tmp.data(), lam/rho, alpha.data());
    // u update
    u += alpha - Dth_tmp;
    // double cur_objective = tf_objective(y, theta, xd, weights, lam, k);

    // Check for convergence
    rr = (Dth_tmp - alpha).norm() / alpha.size();
    ss = rho * Dktv(alpha - alpha_old, k, xd).norm() / dk_mat_sq.rows();
    alpha_old = alpha;
    if (rr < tol && ss < tol) break;
  }
  obj_val = tf_objective(y, theta, xd, weights, lam, k);
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
    int linear_solver = 2) {

  int n = x.size();

  if (lambda[0] < tol / 100 && lambda_max <= 0) {
    lambda_max = get_lambda_max(x, y, weights, k);
  }
  get_lambda_seq(lambda, lambda_max, lambda_min, lambda_min_ratio, nlambda);


  Eigen::MatrixXd theta(n, nlambda);
  Rcpp::NumericVector objective_val(nlambda);
  Rcpp::IntegerVector iters(nlambda);
  Rcpp::IntegerVector dof(nlambda);

  // Use DP solution for k=0.
  if (k == 0) {
    for (int i = 0; i < nlambda; i++) {
      tf_dp_weight(n, y.data(), weights.data(), lambda[i],
                   theta.col(i).data());
      objective_val[i] = tf_objective(y, theta.col(i), x, weights, lambda[i], k);
      dof[i] = calc_degrees_of_freedom(theta.col(i), k);
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
      theta.col(i), alpha.col(i), u, 
      iters[i], objective_val[i], dk_mat,
      dk_mat_sq, lambda[i], max_iter, lambda[i]*rho_scale,
      tol, linear_solver);
    dof[i] = calc_degrees_of_freedom(alpha.col(i), k);
    if (i + 1 < nlambda) {
      theta.col(i + 1) = theta.col(i);
      alpha.col(i + 1) = alpha.col(i);
    }
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("theta") = theta,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("lambda") = lambda,
    Rcpp::Named("tf_objective") = objective_val,
    Rcpp::Named("iters") = iters,
    Rcpp::Named("dof") = dof
  );
  return out;
}


Rcpp::List admm_single_lambda_with_tracking(NumericVector x,
    Eigen::VectorXd& y, const Eigen::ArrayXd& weights, int k,
    double lam, int max_iter, double rho,
    int linear_solver = 2) {

  int n = x.size();

  // Initialize difference matrices and other helper objects
  SparseMatrix<double> penalty_mat = get_penalty_mat(k+1, x);
  SparseMatrix<double> dk_mat = get_dk_mat(k, x, false);
  SparseMatrix<double> dk_mat_sq = dk_mat.transpose() * dk_mat;
  VectorXd wy = (y.array()*weights).matrix();

  // Initialize theta, alpha, u
  VectorXd theta = project_polynomials(x, y, weights, k);
  VectorXd alpha = dk_mat * theta;
  VectorXd u = init_u((theta-y)/rho, x, k, weights);
  VectorXd tmp(n-k);

  if ((linear_solver == 0) & (k != 1)) {
    throw std::invalid_argument("`tridiag` can only be used with k=1.");
  }
  LinearSystem linear_system;
  linear_system.construct(y, weights, k, rho, dk_mat, dk_mat_sq, linear_solver);
  linear_system.compute(linear_solver);

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
        alpha + u, k, x, rho, linear_solver);
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

void LinearSystem::kf_init(const Eigen::VectorXd& y, const Eigen::ArrayXd& weights, int k, double rho, 
              const Eigen::SparseMatrix<double>& dk_mat) {
  int n = y.size();
  wsqrt = weights.sqrt();
  yw = y.array() * wsqrt.array();
  // construct transition matrix A:
  T = MatrixXd::Zero(k, k);
  T.block(1, 0, k - 1, k - 1).diagonal() = VectorXd::Ones(k - 1);
  // construct dense Dseq with only nonzero entries from sparse Dseq:
  Dseq = smat_to_mat2(dk_mat, k);
  // ideally, construct dense Dseq directly:
  //MatrixXd Dseq = b_mat_triplet(k, x, VectorXi::LinSpaced(n - k, 0, n - k - 1));
  // other initialization: 
  RQR = VectorXd::Zero(n - k);
  a1 = VectorXd::Zero(k);
  at = MatrixXd::Zero(k, n + 1);
  P1 = MatrixXd::Zero(k, k);
  Pt = MatrixXd::Zero(k * k, n + 1);
  P1inf = MatrixXd::Identity(k, k);
  Pinf = MatrixXd::Zero(k * k, n + 1);
  Pinf.col(0) = Map<Eigen::VectorXd>(P1inf.data(), k * k);
  // forward: 
  d = 0;
  rankp = k;
  vt_b = 0.0;
  Ft_b = 0.0;
  Finf_b = 0.0;
  vt = VectorXd::Zero(n);
  Ft = VectorXd::Zero(n);
  Finf = VectorXd::Zero(n);
  Kt = MatrixXd::Zero(k, n);
  Kinf = MatrixXd::Zero(k, n);
  Kt_b = VectorXd::Zero(k);
  Kinf_b = VectorXd::Zero(k);
  // re-construct Dseq in which each row saves values for T.row(0) per iterate: 
  s_seq = VectorXd::Ones(n - k);
  s_seq = Dseq.block(0, k, n - k, 1);
  Dseq.conservativeResize(n - k, k);
  MatrixXd firstRow(1, k);
  for (int i = 0; i < n - k; i++) {
    firstRow = -Dseq.row(i) / s_seq(i);
    std::reverse(firstRow.data(), firstRow.data() + k);
    Dseq.row(i) = firstRow;
  }
  RQR = s_seq.array().square().inverse();
  RQR *= 1 / rho; 
  T.row(0) = Dseq.row(0);
  // backward: 
  r = VectorXd::Zero(k);
  r1 = VectorXd::Zero(k);
  rtmp = VectorXd::Zero(k);
  L0 = MatrixXd::Zero(k, k);
  L1 = MatrixXd::Zero(k, k);
  Ptemp = MatrixXd::Zero(k, k);
}

void LinearSystem::kf_iter(const Eigen::VectorXd& c) {
  int n = yw.size();
  int k = a1.size();
  sol = VectorXd::Zero(n);
  // initialize or reset iterator
  d = 0;
  rankp = k;
  // initialize or reset states and covariance matrices
  a1 = VectorXd::Zero(k);
  P1 = MatrixXd::Zero(k, k);
  P1inf = MatrixXd::Identity(k, k);
  Pinf.col(0) = Map<Eigen::VectorXd>(P1inf.data(), k * k);
  while (rankp > 0 && d < n) {
    df1step(yw(d), wsqrt(d), 1, T, RQR(0), a1, P1, P1inf, rankp, vt_b, Ft_b,
            Finf_b, Kt_b, Kinf_b);
    at.col(d + 1) = a1;
    vt(d) = vt_b;
    Ft(d) = Ft_b;
    Finf(d) = Finf_b;
    Kt.col(d) = Kt_b;
    Kinf.col(d) = Kinf_b;
    Pt.col(d + 1) = Map<VectorXd>(P1.data(), k * k);
    Pinf.col(d + 1) = Map<VectorXd>(P1inf.data(), k * k);
    d++;
  }
  for (int i = d; i < n; i++) {
    f1step(yw(i), c(i - d) / s_seq(i - d), wsqrt(i), 1, T, RQR(i - d), a1, P1,
           vt_b, Ft_b, Kt_b);
    vt(i) = vt_b;
    Ft(i) = Ft_b;
    Kt.col(i) = Kt_b;
    at.col(i + 1) = a1;
    Pt.col(i + 1) = Map<VectorXd>(P1.data(), k * k);
    // update transition matrix T for next iterate:
    if (i < n - 1)
      T.row(0) = Dseq.row(i - d + 1);
  }
  for (int i = n - 1; i >= d; i--) {
    P1 = Map<MatrixXd>(Pt.col(i + 1).data(), k, k);
    sol(i) = at.col(i + 1)(0) - P1.row(0) * r;
    L0.col(0) = Kt.col(i) * wsqrt(i) / Ft(i);
    r1 = r - L0.transpose() * r;
    r1(0) -= wsqrt(i) * vt(i) / Ft(i);
    r = T.transpose() * r1;
    // update transition matrix T for next iterate:
    if (i > d)
      T.row(0) = Dseq.row(i - d - 1);
  }
  for (int i = d - 1; i >= 0; i--) {
    P1 = Map<MatrixXd>(Pt.col(i + 1).data(), k, k);
    P1inf = Map<MatrixXd>(Pinf.col(i + 1).data(), k, k);
    a1 = at.col(i + 1);
    sol(i) = a1(0) - P1.row(0) * r - P1inf.row(0) * r1;
    if (i > 0) {
      if (Finf(i) > 0) {
        // simple version w/o mat multiplication:
        L0.col(0) = Kinf.col(i) * wsqrt(i) / Finf(i);
        L1.col(0) = L0.col(0) * Ft(i) / Finf(i);
        L1.col(0) -= Kt.col(i) * wsqrt(i) / Finf(i);
        rtmp = r1 - L0.transpose() * r1 + L1.transpose() * r;
        rtmp(0) -= wsqrt(i) * vt(i) / Finf(i);
        r1 = T.transpose() * rtmp;
        rtmp = r - L0.transpose() * r;
        r = T.transpose() * rtmp;
      } else {
        L1.col(0) = Kt.col(i) * wsqrt(i) / Ft(i);
        rtmp = r - L1.transpose() * r;
        rtmp(0) -= wsqrt(i) * vt(i) / Ft(i);
        r = T.transpose() * rtmp;
        rtmp = r1 - L1.transpose() * r1;
        r1 = T.transpose() * rtmp;
      }
    }
  }
}
