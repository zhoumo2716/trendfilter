// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Numer;

// ---- Define the objective functor ----
class ZObjective : public MFuncGrad {
private:
  const Eigen::VectorXd &W;
  const Eigen::SparseMatrix<double> &dk_mat;
  const Eigen::VectorXd &alpha;
  const Eigen::VectorXd &u;
  double rho;
  int n;
  mutable int eval_count = 0;

public:
  ZObjective(const Eigen::VectorXd &W_,
             const Eigen::SparseMatrix<double> &dk_,
             const Eigen::VectorXd &alpha_,
             const Eigen::VectorXd &u_,
             double rho_, int n_)
    : W(W_), dk_mat(dk_), alpha(alpha_), u(u_), rho(rho_), n(n_) {}

  // --- Compute f(z) and its gradient ---
  double f_grad(Constvec& zcur, Refvec grad) {
    eval_count++;

    // Clip z to prevent overflow in exp()
    Eigen::VectorXd zc = zcur.cwiseMax(-20.0).cwiseMin(20.0);
    Eigen::VectorXd r = dk_mat * zc - alpha - u;

    // Objective: f(z)
    double fval = (W.array() * zc.array().exp()).sum()
      - zc.segment(1, n).sum()
      + 0.5 * rho * r.squaredNorm();

      // Gradient components
      Eigen::VectorXd grad_lik = (W.array() * zc.array().exp()).matrix(); // likelihood term
      grad_lik.segment(1, n).array() -= 1.0;

      Eigen::VectorXd grad_pen = rho * dk_mat.transpose() * r;            // penalty term
      grad = grad_lik + grad_pen;

      // Debug info every 10 evaluations
      if (eval_count % 10 == 0) {
        Rcpp::Rcout << "[ZObjective] eval=" << eval_count
                    << "  f=" << fval
                    << "  ||grad_lik||=" << grad_lik.norm()
                    << "  ||grad_pen||=" << grad_pen.norm()
                    << "  ratio(q/f)=" << grad_pen.norm() / (grad_lik.norm() + 1e-12)
                    << "  ||res||=" << r.norm()
                    << "  z_range=[" << zc.minCoeff() << ", " << zc.maxCoeff() << "]"
                    << std::endl;
        Rcpp::Rcout.flush();
      }

      return fval;
  }
};

// ---- Main optimizer wrapper ----
// [[Rcpp::export]]
Eigen::VectorXd lbfgs_update(
    const Eigen::VectorXd &z_init,
    const Eigen::VectorXd &W,
    int n,
    const Eigen::SparseMatrix<double>& dk_mat,
    const Eigen::VectorXd &alpha,
    const Eigen::VectorXd &u,
    double rho,
    int max_iters = 200,
    double tol = 1e-6)
{
  Eigen::VectorXd z = z_init;
  ZObjective f(W, dk_mat, alpha, u, rho, n);
  double fopt;

  // Run optimizer
  int status = optim_lbfgs(f, z, fopt, max_iters, tol, tol);

  Rcpp::Rcout << "[lbfgs_update] status=" << status
              << "  fopt=" << fopt
              << "  tol=" << tol
              << "  (0 means success)" << std::endl;
  Rcpp::Rcout.flush();

  return z;
}
