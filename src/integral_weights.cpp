#include "integral_weights.h"
#include <RcppEigen.h>
#include <stdexcept>

// [[Rcpp::depends(RcppEigen)]]

// local 3-point Simpson weights on [x_{j-2},x_j]
static void add_panel_weights(
    Eigen::VectorXd &W,
    double x_left, double x_mid, double x_right,
    int idx_left, int idx_mid, int idx_right
) {
  double h1 = x_mid - x_left;
  double h2 = x_right - x_mid;

  double wL = h1/3.0 + h2/6.0 - (h2*h2)/(6.0*h1);
  double wM = (h1*h1)/(6.0*h2) + h1/2.0 + h2/2.0 + (h2*h2)/(6.0*h1);
  double wR = -(h1*h1)/(6.0*h2) + h1/6.0 + h2/3.0;

  W[idx_left]  += wL;
  W[idx_mid]   += wM;
  W[idx_right] += wR;
}

// [[Rcpp::export]]
Eigen::VectorXd compute_integral_weights(
    const Eigen::VectorXd &x,
    double A,
    double B,
    int n
) {
  int dim = n + 2; // x0=A, ..., xn, x_{n+1}=B
  Eigen::VectorXd W = Eigen::VectorXd::Zero(dim);

  // Build extended node vector including endpoints
  Eigen::VectorXd nodes(dim);
  nodes[0] = A;
  nodes.segment(1, n) = x;
  nodes[n+1] = B;

  // Add Simpson panels
  if (n % 2 == 1) {
    // odd n: last panel is (x_{n-1}, x_n, x_{n+1})
    for (int r = 2; r <= n+1; r += 2) {
      add_panel_weights(W, nodes[r-2], nodes[r-1], nodes[r], r-2, r-1, r);
    }
  } else {
    // even n: last panel is (x_{n-2},x_{n-1},x_n) + trapezoid (x_n,x_{n+1})
    for (int r = 2; r <= n; r += 2) {
      add_panel_weights(W, nodes[r-2], nodes[r-1], nodes[r], r-2, r-1, r);
    }
    double h = nodes[n+1] - nodes[n];
    W[n]   += 0.5 * h;
    W[n+1] += 0.5 * h;
  }

  return W;
}
