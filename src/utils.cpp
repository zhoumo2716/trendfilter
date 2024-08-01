#include <stdexcept>
#include <Eigen/Sparse>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <dspline.h>
#include "utils.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(dspline)]]

typedef Eigen::COLAMDOrdering<int> Ord;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* General utilities */
Eigen::SparseMatrix<double> row_scale(
    Eigen::SparseMatrix<double> A,
    Eigen::ArrayXd v) {
  return v.matrix().asDiagonal() * A;
}
Eigen::MatrixXd row_scale(Eigen::MatrixXd A, Eigen::ArrayXd v) {
  return v.matrix().asDiagonal() * A;
}
Eigen::SparseMatrix<double> col_scale(
    Eigen::SparseMatrix<double> A,
    Eigen::ArrayXd v) {
  return A * v.matrix().asDiagonal();
}
Eigen::MatrixXd col_scale(Eigen::MatrixXd A, Eigen::ArrayXd v) {
  return A * v.matrix().asDiagonal();
}

/* Matrix construction */
Eigen::SparseMatrix<double> identity(int n) {
  SparseMatrix<double> Id(n, n);
  Id.setIdentity();
  return Id;
}

Eigen::SparseMatrix<double> diagonal(Eigen::ArrayXd diag) {
  int n = diag.size();
  SparseMatrix<double> D(n, n);
  D.diagonal() = diag;
  return D;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> get_dk_mat(
    int k, NumericVector xd,
    bool tf_weighting) {
  int n = xd.size();
  return dspline::rcpp_b_mat(k, xd, tf_weighting, Rcpp::seq(0, n - k - 1), true);
}

Eigen::SparseMatrix<double> get_penalty_mat(int k, NumericVector xd) {
  int n = xd.size();
  return dspline::rcpp_b_mat(k, xd, true, Rcpp::seq(0, n - k - 1), true);
}

/* Polynomial subspace projection */
Eigen::VectorXd legendre_polynomial(
    Eigen::VectorXd x,
    int k,
    double a,
    double b) {
  ArrayXd xa = 2 * (x.array() - a) / (b - a) - 1;
  if (k == 0) {
    return VectorXd::Ones(x.size());
  } else if (k == 1) {
    return xa.matrix();
  } else if (k == 2) {
    return (1.5*xa.pow(2) - 0.5).matrix();
  } else if (k == 3) {
    return (2.5*xa.pow(3) - 1.5*xa).matrix();
  } else {
    throw std::invalid_argument("`k` must be 0, 1, 2, or 3.");
  }
}

Eigen::MatrixXd polynomial_basis(
    const Eigen::VectorXd& x,
    int k,
    double a = 0.0,
    double b = 1.0) {
  int n = x.size();
  MatrixXd basis_mat(n, k + 1);
  for (int j = 0; j < k + 1; j++) {
    basis_mat.col(j) = legendre_polynomial(x, j, a, b);
  }
  return basis_mat;
}

Eigen::VectorXd project_polynomials(
    const NumericVector& x, // should this be Eigen?
    const VectorXd& y,
    const ArrayXd& weights,
    int k) {
  Eigen::ColPivHouseholderQR<MatrixXd> qr;
  ArrayXd sqrt_weights = weights.sqrt();
  VectorXd x_vec(Rcpp::as<Eigen::VectorXd>(x));
  MatrixXd basis_mat = polynomial_basis(
    x_vec, k, x_vec.minCoeff(), x_vec.maxCoeff()
  );
  // If this isn't accurate enough, can also use SVD.
  qr.compute(row_scale(basis_mat, sqrt_weights));
  VectorXd beta = qr.solve((y.array()*sqrt_weights).matrix());
  VectorXd projection = basis_mat*beta;
  // if (qr.info() > 0) {
  //  std::cerr << "Eigen QR solve returned nonzero exit status.\n";
  // }
  return projection;
}

/* Tridiagonal matrix solve */
Eigen::VectorXd tridiag_forward(
    const Eigen::VectorXd& a,
    const Eigen::VectorXd& b,
    const Eigen::VectorXd& c) {
  int n = a.size();
  Eigen::VectorXd cp(n - 1);

  // Forward sweep part 1
  cp[0] = c[0] / b[0];
  for (int i = 1; i < n - 1; i++) {
    cp[i] = c[i] / (b[i] - a[i]*cp[i - 1]);
  }
  return cp;
}

// Technically, construction of dp is also part of the forward sweep,
// but it makes sense to include it in "backsolve" since d can change
// over iterations.
Eigen::VectorXd tridiag_backsolve(
    const Eigen::VectorXd& a,
    const::VectorXd& b,
    const Eigen::VectorXd& cp,
    const Eigen::VectorXd& d) {
  int n = d.size();
  Eigen::VectorXd dp(n);
  Eigen::VectorXd x(n);

  // Forward sweep part 2
  dp[0] = d[0] / b[0];
  for (int i = 1; i < n; i++) {
    dp[i] = (d[i] - a[i]*dp[i - 1]) / (b[i] - a[i]*cp[i - 1]);
  }

  // Backsolve
  x[n - 1] = dp[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    x[i] = dp[i] - cp[i]*x[i + 1];
  }
  return x;
}

std::tuple<Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd> extract_tridiag(
    Eigen::SparseMatrix<double> A) {
  int n = A.cols();
  VectorXd a(n);
  VectorXd b(n);
  VectorXd c(n);
  // Extract diagonal into b
  b = A.diagonal();
  // Extract (-1)-diagonal into a
  a[0] = 0;
  for (int i = 1; i < n; i++) {
    a[i] = A.coeff(i, i - 1);
  }
  // Extract (+1)-diagonal into c
  for (int i = 0; i < n-1; i++) {
    c[i] = A.coeff(i, i + 1);
  }
  c[n - 1] = 0;
  return std::make_tuple(a, b, c);
}

/* Miscellaneous */
// [[Rcpp::export]]
double get_lambda_max(
    const NumericVector& x,
    const Eigen::VectorXd& y,
    const Eigen::ArrayXd& weights,
    int k) {
  ArrayXd sqrt_weights = weights.sqrt();
  SparseMatrix<double> ck1_mat = get_dk_mat(k + 1, x, true);
  SparseQR<SparseMatrix<double>, Ord> qr;
  qr.compute(col_scale(ck1_mat, sqrt_weights.inverse()).transpose());
  VectorXd u_infty = qr.solve((y.array() * sqrt_weights).matrix());
  return u_infty.lpNorm<Eigen::Infinity>();
}

// [[Rcpp::export]]
int calc_degrees_of_freedom(Eigen::VectorXd const &v, int k, double tol) {
  int dof = k + 1;
  for (int i = 1; i < v.size(); i++) {
    if (abs(v(i) - v(i - 1)) > tol) {
      dof++;
    }
  }
  return dof;
}

// [[Rcpp::export]]
Eigen::VectorXd get_lambda_seq_r(
    Eigen::VectorXd lambda,
    double lambda_max,
    double lambda_min,
    double lambda_min_ratio,
    int n_lambda) {

  get_lambda_seq(lambda, lambda_max, lambda_min, lambda_min_ratio, n_lambda);
  return lambda;
}

void get_lambda_seq(
    Eigen::VectorXd& lambda,
    double lambda_max,
    double lambda_min = -1.0,
    double lambda_min_ratio = 1e-5,
    int n_lambda = 50) {

  if (!(lambda.array() < 1e-12).all()) {
    lambda_min = lambda.minCoeff();
    lambda_max = lambda.maxCoeff();
    n_lambda = lambda.size();
  } else {
    double lmpad = lambda_min_ratio * lambda_max;
    lambda_min = (lambda_min < 0) ? lmpad : lambda_min;
    double ns = static_cast<double>(n_lambda) - 1;
    double p = 0.0;
    lambda[0] = lambda_max;
    if (lambda_min > 1e-20) {
      p = pow(lambda_min / lambda_max, 1 / ns);
      for (int i = 1; i < n_lambda; i++) lambda[i] = lambda[i - 1] * p;
    } else {
      ns -= 1;
      p = pow(lmpad / lambda_max, 1 / ns);
      for (int i = 1; i < n_lambda - 1; i++) lambda[i] = lambda[i - 1] * p;
      lambda[n_lambda - 1] = lambda_min;
    }
  }
}


// ---- Workarounds to access NumericVector utilities in dspline
// ---- Do we alter the function signatures in dspline?

// [[Rcpp::export]]
Eigen::VectorXd Dkv(Eigen::VectorXd v, int k, const NumericVector& xd,
                    bool tf_weighting) {
  Rcpp::NumericVector nv(Rcpp::wrap(v));
  Rcpp::NumericVector out = dspline::rcpp_d_mat_mult(nv, k, xd, tf_weighting, false);
  return Rcpp::as<Eigen::Map<VectorXd> >(out);
}


Eigen::VectorXd Dktv(Eigen::VectorXd v, int k, const NumericVector& xd) {
  Rcpp::NumericVector nv(Rcpp::wrap(v));
  Rcpp::NumericVector out = dspline::rcpp_d_mat_mult(nv, k, xd, false, true);
  return Rcpp::as<Eigen::Map<VectorXd> >(out);
}

// in-place computation of Ptemp given P and A
Eigen::MatrixXd computePtemp(Eigen::MatrixXd A, Eigen::MatrixXd P) {
  int k = P.rows();
  Eigen::MatrixXd temp = A.row(0) * P;
  Eigen::MatrixXd var = A.row(0) * temp.transpose();

  // in-place block replacement starting from the last component:
  P.block(1, 1, k - 1, k - 1).reverse() = P.block(0, 0, k - 1, k - 1).reverse();
  P(0, 0) = var(0, 0);
  temp.conservativeResize(1, k - 1);  // drop the last component
  P.block(0, 1, 1, k - 1) = temp;
  P.block(1, 0, k - 1, 1) = temp.transpose();
  return P;
}

// [[Rcpp::export]]
Eigen::MatrixXd smat_to_mat(const Eigen::SparseMatrix<double>& sparseMat, int k, bool equal_spaced) {
  int rows = sparseMat.rows(); // n-k
  Eigen::MatrixXd denseMat(rows, k + 1);
  std::vector<double> rowNonzeros;
  // Iterate over nonzero coefficients in each row of the sparse matrix
  for (int i = 0; i < sparseMat.outerSize(); ++i) {
    std::vector<double> rowNonzeros;
    for (SparseMatrix<double>::InnerIterator it(sparseMat, i); it; ++it) 
      rowNonzeros.push_back(it.value());
    int m = rowNonzeros.size();
    for (int j = 0; j < m; j++) {
      if (i < rows) 
        denseMat(i - j, j) = rowNonzeros[m - 1 - j];
      else  
        denseMat(rows - 1 - j, j + i - rows + 1) = rowNonzeros[m - 1 - j];
    }
    if (equal_spaced && i == k) {
      denseMat.conservativeResize(1, k + 1);
      return denseMat;
    } 
  }
  return denseMat;
}

void f1step(double y, double c, double Z, double H, const Eigen::MatrixXd& A, 
  double RQR, Eigen::VectorXd& a, Eigen::MatrixXd& P, double& vt, double& Ft,
  Eigen::VectorXd& Kt) {
  VectorXd a_temp = A * a;
  a_temp(0) += c;
  MatrixXd Ptemp = computePtemp(A, P);  // A * P * A.transpose();
  Ptemp(0, 0) += RQR;

  vt = y - Z * a_temp(0);
  Ft = pow(Z, 2) * Ptemp(0, 0) + H;
  Kt = Ptemp.col(0) * Z;
  a = a_temp + Kt * vt / Ft;
  P = Ptemp - Kt * Z * Ptemp.row(0) / Ft;

  // symmetrize
  Ptemp = P;
  Ptemp += P.transpose();
  P = Ptemp / 2;

  // Some entries of P can be _really_ small in magnitude, set them to 0
  // but maintain symmetric / posdef
  for (int i = 0; i < P.rows(); i++) {
    for (int j = 0; j <= i; j++) {
      if (abs(P(i, j)) < 1e-30) {
        if (i == j) {
          P.row(i).setZero();
          P.col(i).setZero();
        } else {
          P(i, j) = 0;
          P(j, i) = 0;
        }
      }
    }
  }
}

void df1step(double y, double Z, double H, const Eigen::MatrixXd& A, double RQR,
  Eigen::VectorXd& a, Eigen::MatrixXd& P, Eigen::MatrixXd& Pinf, int& rankp,
  double& vt, double& Ft, double& Finf, Eigen::VectorXd& Kt, Eigen::VectorXd& Kinf) {
  double tol = Eigen::NumTraits<double>::epsilon();
  tol = std::sqrt(tol);
  int k = a.size();
  MatrixXd Ptemp(k, k);

  a = A * a;
  P = computePtemp(A, P);  // A * P * A.transpose();
  P(0, 0) += RQR;
  Pinf = computePtemp(A, Pinf);  // A * Pinf * A.transpose();

  vt = y - Z * a(0);
  Kt = P.col(0) * Z;
  Ft = Z * Kt(0) + H;
  Kinf = Pinf.col(0) * Z;
  Finf = Z * Kinf(0);

  if (Finf > tol) {  // should always happen
    a += vt * Kinf / Finf;
    P += Ft * Kinf * Kinf.transpose() / pow(Finf, 2);
    P -= Kt * Kinf.transpose() / Finf + Kinf * Kt.transpose() / Finf;
    Pinf -= Kinf * Kinf.transpose() / Finf;
    rankp--;
  } else {  // should never happen
    Finf = 0;
    if (Ft > tol) {
      a += vt * Kt / Ft;
      P -= Kt * Kt.transpose() / Ft;
    }
  }
  if (Ft < tol)
    Ft = 0;

  // symmetrize
  Ptemp = P;
  Ptemp += P.transpose();
  P = Ptemp / 2;
  Ptemp = Pinf;
  Ptemp += Pinf.transpose();
  Pinf = Ptemp / 2;

  // Fix possible negative definiteness, should never happen
  for (int i = 0; i < k; i++) {
    if (P(i, i) < 0) {
      P.row(i).setZero();
      P.col(i).setZero();
    }
    if (Pinf(i, i) < 0) {
      Pinf.row(i).setZero();
      Pinf.col(i).setZero();
    }
  }
}
