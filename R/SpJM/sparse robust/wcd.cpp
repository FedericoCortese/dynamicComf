#include <Rcpp.h>
using namespace Rcpp;

// helper to compute median of a std::vector<double>
// includes the diagonal zero entry as in your R code
double median_vec(std::vector<double>& v) {
  std::sort(v.begin(), v.end());
  int n = v.size();
  if (n % 2 == 1) {
    return v[(n - 1) / 2];
  } else {
    return 0.5 * (v[n/2 - 1] + v[n/2]);
  }
}

// [[Rcpp::export]]
NumericMatrix WCD(const IntegerVector& s,
                      const NumericMatrix& Y,
                      int K) {
  int TT = Y.nrow();
  int P  = Y.ncol();
  
  // output: K x P
  NumericMatrix wcd(K, P);
  
  // compute sk[p] = IQR(Y[,p]) / 1.35
  std::vector<double> sk(P);
  for (int p = 0; p < P; ++p) {
    std::vector<double> col(TT);
    for (int i = 0; i < TT; ++i) col[i] = Y(i,p);
    std::sort(col.begin(), col.end());
    // type‐1 quantiles at 25% and 75%
    double q1 = col[(int)std::floor((TT - 1) * 0.25)];
    double q3 = col[(int)std::floor((TT - 1) * 0.75)];
    sk[p] = (q3 - q1) / 1.35;
    if (sk[p] == 0.0) sk[p] = 1.0;
  }
  
  // for each cluster i = 1..K
  for (int ci = 1; ci <= K; ++ci) {
    // collect row‐indices belonging to cluster ci
    std::vector<int> rows;
    for (int i = 0; i < TT; ++i) {
      if (s[i] == ci) rows.push_back(i);
    }
    int n = rows.size();
    if (n < 2) {
      // leave wcd(ci-1, *) as zeros
      continue;
    }
    
    // for each dimension p
    for (int p = 0; p < P; ++p) {
      // build the within‐cluster distance matrix for dimension p
      NumericMatrix mat(n, n);
      for (int a = 0; a < n; ++a) {
        mat(a,a) = 0.0;
        for (int b = a+1; b < n; ++b) {
          double d = std::abs(Y(rows[a],p) - Y(rows[b],p)) / sk[p];
          mat(a,b) = mat(b,a) = d;
        }
      }
      // compute row‐medians (including the zero on the diagonal)
      double sum_meds = 0.0;
      std::vector<double> tmp(n);
      for (int a = 0; a < n; ++a) {
        for (int j = 0; j < n; ++j) tmp[j] = mat(a,j);
        sum_meds += median_vec(tmp);
      }
      // average of the medians
      wcd(ci-1, p) = sum_meds / n;
    }
  }
  
  return wcd;
}
