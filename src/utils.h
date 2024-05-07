#ifndef UTILS_H
#define UTILS_H

using Rcpp::NumericVector;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

/* Matrix construction */
Eigen::SparseMatrix<double> identity(int n);
Eigen::SparseMatrix<double> diagonal(Eigen::ArrayXd diag);
Eigen::SparseMatrix<double> get_dk_mat(int k, NumericVector xd,
    bool tf_weighting);
Eigen::SparseMatrix<double> get_penalty_mat(int k, NumericVector xd);
Eigen::VectorXd legendre_polynomial(Eigen::VectorXd, int k, double a,
    double b);

/* Polynomial subspace projection */
Eigen::VectorXd project_polynomials(const NumericVector& x, const VectorXd& y,
    const Eigen::ArrayXd& weights, int k);

/* Tridiagonal matrix solve */
std::tuple<Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd> extract_tridiag(
    Eigen::SparseMatrix<double> A);
Eigen::VectorXd tridiag_forward(const Eigen::VectorXd& a,
        const Eigen::VectorXd& b, const Eigen::VectorXd& c);
Eigen::VectorXd tridiag_backsolve(
        const Eigen::VectorXd& a, const::VectorXd& b,
            const Eigen::VectorXd& cp, const Eigen::VectorXd& d);

/* Lambda sequence */
double get_lambda_max(const NumericVector& x, const Eigen::VectorXd& y,
                      const Eigen::ArrayXd& sqrt_weights, int k);
void get_lambda_seq(Eigen::VectorXd& lambda, double lambda_max,
                    double lambda_min, double lambda_min_ratio, int n_lambda);
Eigen::VectorXd get_lambda_seq_r(Eigen::VectorXd lambda, double lambda_max,
                                 double lambda_min, double lambda_min_ratio,
                                 int n_lambda);


#endif

