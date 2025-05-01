#ifndef SOFTTHRESHOLD_H
#define SOFTTHRESHOLD_H

#include <RcppEigen.h>

// Declaration of the soft-thresholding function.
// This function applies the soft-thresholding operator elementwise to a numeric vector.
// For each component: S_gamma(x)_i = sign(x_i) * max(|x_i| - gamma, 0)
Eigen::VectorXd softThreshold(const Eigen::VectorXd &x, double gamma);

#endif
