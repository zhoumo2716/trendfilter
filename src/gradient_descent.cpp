#include <cmath>
#include <Rcpp.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "gradient_descent.h"

using Eigen::VectorXd;
using Eigen::SparseMatrix;

// Gradient-descent z-update with Armijo backtracking (no Hessian)
// Objective: phi(z) = sum_i W_i * exp(z_i) - sum_{i=1}^n z_i + (rho/2)*||Dk*z - alpha - u||^2
VectorXd gd_update(
    const VectorXd &z_init,
    const VectorXd &W,
    int n,
    const SparseMatrix<double>& dk_mat,
    //const SparseMatrix<double>& dk_mat_sq_unused,
    const VectorXd &alpha,
    const VectorXd &u,
    double rho,
    int max_iters,
    double tol)
{
  const int dim = n + 2;
  VectorXd z = z_init;

  // --- keep exp() safe ---
  const double z_lo = -40.0;   // exp(-40) ~ 4e-18
  const double z_hi =  40.0;   // exp(40)  ~ 2e17
  // Objective
  auto phi = [&](const VectorXd& zcur) -> double {
    // term1 = sum W_i * exp(z_i)   (use clipped z to avoid Inf/NaN)
    VectorXd zc = zcur.cwiseMax(z_lo).cwiseMin(z_hi);
    double term1 = (W.array() * zc.array().exp()).sum();
    // term2 = sum_{i=1..n} z_i
    double term2 = zcur.segment(1, n).sum();
    // term3 = (rho/2) * ||Dk z - alpha - u||^2
    VectorXd r = dk_mat * zcur - alpha - u;
    double term3 = 0.5 * rho * r.squaredNorm();
    return term1 - term2 + term3;
  };
  // Backtracking parameters
  const double c1    = 1e-4;    // Armijo slope fraction
  const double beta  = 0.5;     // step shrink
  const double t_min = 1e-12;   // min step to avoid stalling
  double fval = phi(z);


  for (int it = 0; it < max_iters; ++it) {
    // ---- Gradient = grad_f + grad_q ----
    // grad_f (likelihood part)
    VectorXd zc = z.cwiseMax(z_lo).cwiseMin(z_hi);
    Eigen::ArrayXd expz = zc.array().exp();
    VectorXd grad_f = (W.array() * expz).matrix();  // size dim
    grad_f.segment(1, n).array() -= 1.0;            // subtract 1 on interior indices

    // grad_q (ADMM quadratic part)
    VectorXd r = dk_mat * z - alpha - u;            // size dk_mat.rows()
    VectorXd grad_q = rho * dk_mat.transpose() * r; // size dim

    // total gradient
    VectorXd grad = grad_f + grad_q;

    if (!grad.allFinite())
      Rcpp::stop("[gd_update] gradient has non-finite entries (NaN/Inf).");

    // First-order optimality
    if (grad.lpNorm<Eigen::Infinity>() < tol)
      break;

    // Descent direction & slope
    VectorXd d = -grad;
    double gTd = grad.dot(d); // = -||grad||^2 <= 0
    if (gTd >= 0) {
      d = -grad;
      gTd = -grad.squaredNorm();
    }
    // ---- Backtracking line search (Armijo) ----
    double t = 1.0;
    bool accepted = false;
    while (t >= t_min) {
      VectorXd z_new = z + t * d;
      if (!z_new.allFinite()) { t *= beta; continue; }
      double f_new = phi(z_new);
      if (std::isfinite(f_new) && f_new <= fval + c1 * t * gTd) {
        z = z_new;
        fval = f_new;
        accepted = true;
        break;
      }
      t *= beta;
    }

    if (!accepted) {
      // Fallback: tiny step; if still not finite, bail.
      VectorXd z_new = z + t_min * d;
      if (!z_new.allFinite())
        Rcpp::stop("[gd_update] z became non-finite after fallback step.");
      z = z_new;
      fval = phi(z);
    }

    // Step-size convergence (cheap early stop)
    if ((d * std::max(1e-12, std::min(1.0, 0.5))).norm() < tol)
      break;
  }

  return z;
}
