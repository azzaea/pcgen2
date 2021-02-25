

// [[Rcpp::export]]
double fcn(Rcpp::NumericVector v){
  auto w = Rcpp::as<vector<double>>(v);
  double y = doSomething(w);
  return y;
}



double doSomething(const std::vector<double>& v){
  double ret = 0.1;
  return ret;
}
