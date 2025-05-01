#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
Eigen::VectorXd softThreshold(const Eigen::VectorXd &x, double gamma) {
  int n = x.size();
  Eigen::VectorXd result(n);
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    if (xi > gamma) {
      result(i) = xi - gamma;
    } else if (xi < -gamma) {
      result(i) = xi + gamma;
      } else {
      result(i) = 0.0;
        }
  }
  return result;
  }
