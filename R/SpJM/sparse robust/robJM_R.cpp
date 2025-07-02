#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <set>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gower_dist(const NumericMatrix& Y,
                         const NumericMatrix& mu) {
  int n = Y.nrow();     // number of observations in Y
  int p = Y.ncol();     // number of variables (columns)
  int m = mu.nrow();    // number of rows in mu
  
  // Compute s_p: the range of each column of Y
  NumericVector s_p(p);
  for(int j = 0; j < p; ++j) {
    double mn = Y(0, j), mx = Y(0, j);
    for(int i = 1; i < n; ++i) {
      if (Y(i, j) < mn) mn = Y(i, j);
      if (Y(i, j) > mx) mx = Y(i, j);
    }
    s_p[j] = mx - mn;
    if (s_p[j] == 0.0) s_p[j] = 1.0;  // avoid division by zero
  }
  
  // Allocate output matrix V (n x m)
  NumericMatrix V(n, m);
  
  // Compute the Gower distance (mean of standardized abs differences)
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      double acc = 0.0;
      for(int k = 0; k < p; ++k) {
        acc += std::abs(Y(i, k) - mu(j, k)) / s_p[k];
      }
      V(i, j) = acc / p;  // divide by number of variables to get the mean
    }
  }
  
  return V;
}


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
NumericVector LOF_cpp(const NumericMatrix& D, int knn) {
  int n = D.nrow();
  if (D.ncol() != n) stop("D must be a square matrix");
  if (knn < 1 || knn >= n) stop("knn must satisfy 1 <= knn < n");
  
  // 1. Compute k-distance and neighbor sets N_k(i)
  std::vector<double> kdist(n);
  std::vector<std::vector<int>> neigh(n);
  for (int i = 0; i < n; ++i) {
    // collect distances to all other points
    std::vector<double> dists;
    dists.reserve(n - 1);
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      dists.push_back(D(i, j));
    }
    std::sort(dists.begin(), dists.end());
    kdist[i] = dists[knn - 1];  // k-th smallest
    // define neighborhood N_k(i)
    for (int j = 0; j < n; ++j) {
      if (j != i && D(i, j) <= kdist[i]) {
        neigh[i].push_back(j);
      }
    }
  }
  
  // 2. Compute local reachability density (lrd)
  std::vector<double> lrd(n);
  for (int i = 0; i < n; ++i) {
    double sum_reach = 0.0;
    for (int j : neigh[i]) {
      double reach_dist = std::max(kdist[j], D(i, j));
      sum_reach += reach_dist;
    }
    if (sum_reach <= 0.0) {
      lrd[i] = R_PosInf;
    } else {
      lrd[i] = neigh[i].size() / sum_reach;
    }
  }
  
  // 3. Compute LOF for each point
  NumericVector lof(n);
  for (int i = 0; i < n; ++i) {
    double sum_ratio = 0.0;
    for (int j : neigh[i]) {
      sum_ratio += lrd[j] / lrd[i];
    }
    lof[i] = sum_ratio / neigh[i].size();
  }
  
  return lof;
}

// [[Rcpp::export]]
NumericMatrix WCD_old(const IntegerVector& s,
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

// [[Rcpp::export]]
NumericMatrix WCD(const IntegerVector& s,
                  const NumericMatrix& Y,
                  int K,
                  Nullable<IntegerVector> feat_type = R_NilValue) {
  int TT = Y.nrow();
  int P  = Y.ncol();
  
  // 1. Handle NULL feat_type -> all continuous
  IntegerVector ft;
  if (feat_type.isNotNull()) {
    ft = feat_type.get();
    if ((int)ft.size() != P)
      stop("feat_type must have length = ncol(Y)");
  } else {
    ft = IntegerVector(P, 0);
  }
  
  // 2. Compute scale for continuous, dummy for others
  std::vector<double> sk(P);
  for (int p = 0; p < P; ++p) {
    if (ft[p] == 0) {
      std::vector<double> col(TT);
      for (int i = 0; i < TT; ++i) col[i] = Y(i, p);
      std::sort(col.begin(), col.end());
      double q1 = col[(int)std::floor((TT - 1) * 0.25)];
      double q3 = col[(int)std::floor((TT - 1) * 0.75)];
      sk[p] = (q3 - q1) / 1.35;
      if (sk[p] == 0.0) sk[p] = 1.0;
    } else {
      sk[p] = 1.0;
    }
  }
  
  // 3. Precompute ordinal ranks and levels
  std::vector<std::vector<int>> ord_ranks(P);
  std::vector<int> M(P, 1);
  for (int p = 0; p < P; ++p) {
    if (ft[p] == 2) {
      std::vector<double> col(TT);
      for (int i = 0; i < TT; ++i) col[i] = Y(i, p);
      std::vector<double> levels = col;
      std::sort(levels.begin(), levels.end());
      levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
      M[p] = std::max((int)levels.size(), 1);
      ord_ranks[p].resize(TT);
      for (int i = 0; i < TT; ++i) {
        int idx = std::lower_bound(levels.begin(), levels.end(), Y(i, p)) - levels.begin();
        ord_ranks[p][i] = idx;
      }
    }
  }
  
  // 4. Prepare output: K x P
  NumericMatrix wcd(K, P);
  
  // 5. For each cluster
  for (int ci = 1; ci <= K; ++ci) {
    // collect rows in cluster ci
    std::vector<int> rows;
    for (int i = 0; i < TT; ++i) if (s[i] == ci) rows.push_back(i);
    int n = rows.size();
    if (n < 2) continue;
    
    // compute per-dimension WCD
    for (int p = 0; p < P; ++p) {
      // build n x n distance matrix for this feature
      NumericMatrix mat(n, n);
      for (int a = 0; a < n; ++a) {
        mat(a,a) = 0.0;
        for (int b = a+1; b < n; ++b) {
          double diff;
          if (ft[p] == 0) {
            diff = std::abs(Y(rows[a], p) - Y(rows[b], p)) / sk[p];
          } else if (ft[p] == 1) {
            diff = (Y(rows[a], p) != Y(rows[b], p)) ? 1.0 : 0.0;
          } else {
            double denom = double(M[p] - 1);
            diff = denom > 0
            ? std::abs(ord_ranks[p][rows[a]] - ord_ranks[p][rows[b]]) / denom
            : 0.0;
          }
          mat(a,b) = mat(b,a) = diff;
        }
      }
      // compute average row-median
      double sum_meds = 0.0;
      std::vector<double> tmp(n);
      for (int a = 0; a < n; ++a) {
        for (int j = 0; j < n; ++j) tmp[j] = mat(a,j);
        sum_meds += median_vec(tmp);
      }
      wcd(ci-1, p) = sum_meds / n;
    }
  }
  
  return wcd;
}


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix weight_inv_exp_dist(const NumericMatrix& Y,
                                  const IntegerVector& s,
                                  const NumericMatrix& W,
                                  double zeta,  // tuning parameter
                                  Nullable<IntegerVector> feat_type = R_NilValue
) {
  int TT = Y.nrow(), P = Y.ncol();
  
  // 1. Handle NULL feat_type -> all continuous
  IntegerVector ft;
  if (feat_type.isNotNull()) {
    ft = feat_type.get();
    if ((int)ft.size() != P)
      stop("feat_type must have length = ncol(Y)");
  } else {
    ft = IntegerVector(P, 0);
  }
  
  // 2. Compute scale for continuous, dummy for others
  std::vector<double> sk(P);
  for (int p = 0; p < P; ++p) {
    if (ft[p] == 0) {
      std::vector<double> col(TT);
      for (int i = 0; i < TT; ++i) col[i] = Y(i, p);
      std::sort(col.begin(), col.end());
      double q1 = col[(int)std::floor((TT - 1) * 0.25)];
      double q3 = col[(int)std::floor((TT - 1) * 0.75)];
      sk[p] = (q3 - q1) / 1.35;
      if (sk[p] == 0.0) sk[p] = 1.0;
    } else {
      sk[p] = 1.0;
    }
  }
  
  // 3. Precompute ordinal ranks and level counts
  std::vector<std::vector<int>> ord_ranks(P);
  std::vector<int> M(P, 1);
  for (int p = 0; p < P; ++p) {
    if (ft[p] == 2) {
      std::vector<double> col(TT);
      for (int i = 0; i < TT; ++i) col[i] = Y(i, p);
      std::vector<double> levels = col;
      std::sort(levels.begin(), levels.end());
      levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
      int mp = levels.size();
      M[p] = mp > 1 ? mp : 1;
      ord_ranks[p].resize(TT);
      for (int i = 0; i < TT; ++i) {
        int idx = std::lower_bound(levels.begin(), levels.end(), Y(i, p)) - levels.begin();
        ord_ranks[p][i] = idx;
      }
    }
  }
  
  // 4. Build distance matrix
  NumericMatrix D(TT, TT);
  for (int i = 0; i < TT; ++i) {
    for (int j = i + 1; j < TT; ++j) {
      double sum_exp = 0.0;
      int si = s[i] - 1, sj = s[j] - 1;
      for (int p = 0; p < P; ++p) {
        double diff;
        if (ft[p] == 0) {
          diff = std::abs(Y(i, p) - Y(j, p)) / sk[p];
        } else if (ft[p] == 1) {
          diff = (Y(i, p) != Y(j, p)) ? 1.0 : 0.0;
        } else {
          double denom = double(M[p] - 1);
          diff = denom > 0
          ? std::abs(ord_ranks[p][i] - ord_ranks[p][j]) / denom
          : 0.0;
        }
        double wmax = std::max(W(si, p), W(sj, p));
        sum_exp += wmax * std::exp(-diff / zeta);
      }
      D(i, j) = D(j, i) = -zeta * std::log(std::max(sum_exp, 1e-16));
    }
  }
  
  return D;
}





// [[Rcpp::export]]
NumericMatrix weight_inv_exp_dist_old(const NumericMatrix& Y,
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