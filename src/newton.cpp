#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>
#include <iostream>

#include "newton.h"

// Perform Newton update for z-subproblem inside ADMM
Eigen::VectorXd newton_update(
    const Eigen::VectorXd &z_init,
    const Eigen::VectorXd &W,
    int n,
    double A, double B,
    const Eigen::MatrixXd &dk_mat,
    const Eigen::SparseMatrix<double>& dk_mat_sq,
    const Eigen::VectorXd &alpha,
    const Eigen::VectorXd &u,
    double rho,
    int max_iters,
    double tol
) {
  int dim = n + 2;
  Eigen::VectorXd z = z_init;

  const double ZMAX = 40.0, ZMIN = -40.0;

  for (int it = 0; it < max_iters; ++it) {
    Eigen::ArrayXd expz = z.array().min(ZMAX).max(ZMIN).exp();

    // --- Likelihood gradient + Hessian ---
    Eigen::VectorXd s = (W.array() * expz).matrix();

    Eigen::VectorXd grad_f = s;
    for (int i = 1; i <= n; i++) grad_f[i] -= 1.0;
    Eigen::VectorXd diagH_f = s;

    // --- ADMM quadratic term ---
    Eigen::VectorXd residual = dk_mat * z - alpha - u;
    Eigen::VectorXd grad_q = rho * dk_mat.transpose() * residual;
    Eigen::SparseMatrix<double> H_q = rho * dk_mat_sq;

    // --- Combine ---
    Eigen::VectorXd grad = grad_f + grad_q;
    Eigen::SparseMatrix<double> H = H_q;
    H.diagonal().array() += diagH_f.array();

    // --- Newton step ---
    Eigen::VectorXd delta_z;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(H);
    bool ok = (solver.info() == Eigen::Success);
    if (ok) {
      delta_z = solver.solve(-grad);
      ok = (solver.info() == Eigen::Success) && delta_z.allFinite();
    }
    if (!ok) {
      std::cerr << "[Newton] Hessian solve failed, fallback to diagonal.\n";
      Eigen::VectorXd Hdiag = H.diagonal();
      const double eps = 1e-6;
      for (int i = 0; i < Hdiag.size(); ++i)
        if (!std::isfinite(Hdiag[i]) || Hdiag[i] < eps) Hdiag[i] = eps;
        delta_z = (-grad).array() / Hdiag.array();
    }

    // --- Check convergence ---
    if (delta_z.norm() < tol) break;
    z += delta_z;
  }

  // ============================================================
  // === Print Final Diagnostic Info (like L-BFGS version) ======
  // ============================================================

  Eigen::ArrayXd expz_final = z.array().min(ZMAX).max(ZMIN).exp();
  Eigen::VectorXd s_final = (W.array() * expz_final).matrix();

  Eigen::VectorXd grad_lik = s_final;
  for (int i = 1; i <= n; i++) grad_lik[i] -= 1.0;

  Eigen::VectorXd residual_final = dk_mat * z - alpha - u;
  Eigen::VectorXd grad_pen = rho * dk_mat.transpose() * residual_final;

  double fval_final = (W.array() * expz_final).sum()
    - z.segment(1, n).sum()
    + 0.5 * rho * residual_final.squaredNorm();

    std::cout << "[Newton final] f=" << fval_final
              << "  ||grad_lik||=" << grad_lik.norm()
              << "  ||grad_pen||=" << grad_pen.norm()
              << "  z_range=[" << z.minCoeff() << ", " << z.maxCoeff() << "]"
              << std::endl;

    // ============================================================
    return z;
}
