#include "NonmemberCppFcns.h"
#include <vector>
#include <Rcpp.h>

// This seems to be optional...at least for Windows/Rtools:
// [[Rcpp::plugins(cpp17)]]

// Nonmember Function Interfaces:

// [[Rcpp::export]]
int rAdd(double x, double y)
{
  // double add(double x, double y)
  // Call the add(.) function in the reusable C++ code base:
  return add(x, y);
}

// [[Rcpp::export]]
Rcpp::NumericVector rSortVec(Rcpp::NumericVector v)
{
  // vector<double> sortVec(vector<double> v)
  // Transfer data from NumericVector to std::vector<double>
  auto stlVec = Rcpp::as<std::vector<double>>(v);

  // Call the reusable sortVec(.) function, with the expected
  // std::vector<double> argument:
  stlVec = sortVec(stlVec);

  // Reassign the results from the vector<double> return object
  // to the same NumericVector v, using Rcpp::wrap(.):
  v = Rcpp::wrap(stlVec);

  // Return as an Rcpp::NumericVector:
  return v;
}

// C++17 example:
// [[Rcpp::export]]
int rProdLcmGcd(int m, int n)
{
  // int prodLcmGcd(int m, int n)
  return prodLcmGcd(m, n);
}


