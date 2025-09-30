#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>
#include <iostream>

// Perform Newton update for z-subproblem inside ADMM
Eigen::VectorXd newton_update(
    const Eigen::VectorXd &z_init,             // starting point (warm start)
    const std::vector<double> &x,              // knot locations x1,...,xn
    const std::vector<double> &W,              // integration weights W_0,...,W_{n+1}
    int n,                                     // number of interior nodes
    double A, double B,                        // integration bounds
    const Eigen::MatrixXd &dk_mat,             // D^(k) matrix
    const Eigen::SparseMatrix<double>& dk_mat_sq, // (D^(k))^T D^(k)
    const Eigen::VectorXd &alpha,              // α^(t)
    const Eigen::VectorXd &u,                  // u^(t)
    double rho,                                // ADMM penalty parameter
    int max_iters = 10,                        // max Newton iterations
    double tol = 1e-6                          // stopping tolerance
) {
    // dimension is n+2: nodes x0=A, ..., xn, x_{n+1}=B
    int dim = n + 2;
    Eigen::VectorXd z = z_init;

    for (int it = 0; it < max_iters; ++it) {
        // --- Likelihood gradient + Hessian ---
        Eigen::ArrayXd expz = z.array().exp();
        Eigen::VectorXd s(dim);
        for (int i = 0; i < dim; i++) {
            s[i] = W[i] * expz[i];
        }

        // gradient of f(z)
        Eigen::VectorXd grad_f = s;
        for (int i = 1; i <= n; i++) {
            grad_f[i] -= 1.0;  // subtract derivative of -sum_{i=1}^n z_i
        }

        // Hessian diagonal (likelihood)
        Eigen::VectorXd diagH_f = s;

        // --- ADMM quadratic term ---
        Eigen::VectorXd residual = dk_mat * z - alpha - u;
        Eigen::VectorXd grad_q = rho * dk_mat.transpose() * residual;
        Eigen::SparseMatrix<double> H_q = rho * dk_mat_sq;

        // --- Full gradient and Hessian ---
        Eigen::VectorXd grad = grad_f + grad_q;

        Eigen::SparseMatrix<double> H(dim, dim);
        H = H_q;  // start with rho D^T D
        for (int i = 0; i < dim; i++) {
            H.coeffRef(i,i) += diagH_f[i]; // add diagonal from likelihood
        }

        // --- Newton step ---
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(H);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Newton solver factorization failed.\n";
            break;
        }

        Eigen::VectorXd delta_z = solver.solve(-grad);

        // check convergence
        if (delta_z.norm() < tol) {
            break;
        }

        // update step (τ = 1, could backtrack if needed)
        z += delta_z;
    }

    return z;
}
