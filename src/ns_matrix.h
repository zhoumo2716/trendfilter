#ifndef NS_MATRIX_H  // Prevents multiple inclusions
#define NS_MATRIX_H

#include <Rcpp.h>  // Include Rcpp dependency

using namespace Rcpp; // Optional (can be removed if using explicit scoping)

// Declare the function(s) you implemented in `ns_matrix.cpp`
NumericMatrix ns_matrix(NumericVector x, int m);

#endif // NS_MATRIX_H
