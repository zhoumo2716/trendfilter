#ifndef INTEGRAL_WEIGHTS_H
#define INTEGRAL_WEIGHTS_H

#include <vector>

// Compute integration weights W[0..n+1] for nodes x0=A, ..., xn, x_{n+1}=B.
// - x: vector of interior knots {x1,...,xn}, length = n
// - A, B: integration bounds
// Returns: vector<double> of length n+2 with weights W0,...,W_{n+1}.
std::vector<double> compute_integral_weights(
    const std::vector<double> &x,
    double A,
    double B,
    int n
);

#endif // INTEGRAL_WEIGHTS_H
