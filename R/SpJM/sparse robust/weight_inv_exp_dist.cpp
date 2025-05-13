#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix weight_inv_exp_dist(const NumericMatrix& Y,
                                      const IntegerVector& s,
                                      const NumericMatrix& W,
                                      double zeta) {
  int TT = Y.nrow();
  int P  = Y.ncol();
  
  // 1. Compute sk = IQR(col)/1.35 (approximate, type-1 quantiles)
  std::vector<double> sk(P);
  for (int p = 0; p < P; ++p) {
    std::vector<double> col(TT);
    for (int i = 0; i < TT; ++i) col[i] = Y(i, p);
    std::sort(col.begin(), col.end());
    // simple quartiles at 25% and 75% positions
    double q1 = col[(int)std::floor((TT - 1) * 0.25)];
    double q3 = col[(int)std::floor((TT - 1) * 0.75)];
    sk[p] = (q3 - q1) / 1.35;
    if (sk[p] == 0.0) sk[p] = 1.0;  // avoid divide-by-zero
  }
  
  // 2. Prepare output matrix
  NumericMatrix D(TT, TT);
  
  // 3. Loop over all pairs i < j
  for (int i = 0; i < TT; ++i) {
    for (int j = i + 1; j < TT; ++j) {
      double sum_exp = 0.0;
      int si = s[i] - 1; // R→C index
      int sj = s[j] - 1;
      for (int p = 0; p < P; ++p) {
        double diff = std::abs(Y(i, p) - Y(j, p)) / sk[p];
        double wmax = std::max(W(si, p), W(sj, p));
        sum_exp += wmax * std::exp(-diff / zeta);
      }
      // avoid log(0)
      double val = -zeta * std::log(std::max(sum_exp, 1e-16));
      D(i, j) = D(j, i) = val;
    }
  }
  
  return D;
}
