// [[Rcpp::depends(RcppArmadillo)]]  // if you need arma, otherwise remove
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


// [[Rcpp::export]]
IntegerVector initialize_states(const NumericMatrix& Y, int K) {
  int n = Y.nrow(), p = Y.ncol();
  // (1) pick the first centroid uniformly
  int init_idx = floor(R::runif(0, n));          // in [0, n-1]
  std::vector<int> cent_idx;
  cent_idx.reserve(K);
  cent_idx.push_back(init_idx);
  
  // (2) compute distances of every point to that first centroid
  //    we can call your gower_dist(Y, single_row_matrix)
  NumericMatrix first_cent(1, p);
  for(int j = 0; j < p; ++j)
    first_cent(0, j) = Y(init_idx, j);
  
  NumericMatrix D1 = gower_dist(Y, first_cent);  // n × 1
  NumericVector closest_dist = D1(_, 0);         // extract as length-n
  
  // (3) sample the remaining K−1 centroids, *w.p.* ∝ these fixed distances
  for(int i = 1; i < K; ++i) {
    double sumd = std::accumulate(closest_dist.begin(), closest_dist.end(), 0.0);
    NumericVector prob = closest_dist / sumd;
    
    // draw u ∼ Uniform(0,1) and pick the index
    double u = R::runif(0, 1);
    double csum = 0;
    int next_idx = n - 1;
    for(int t = 0; t < n; ++t) {
      csum += prob[t];
      if (u <= csum) { next_idx = t; break; }
    }
    cent_idx.push_back(next_idx);
    // **note**: we do *not* update closest_dist inside the loop, 
    // exactly as in your R code.
  }
  
  // (4) build the centroid matrix of size K×p
  NumericMatrix centroids(K, p);
  for(int i = 0; i < K; ++i)
    for(int j = 0; j < p; ++j)
      centroids(i, j) = Y(cent_idx[i], j);
  
  // (5) assign each observation to its nearest centroid
  NumericMatrix Dk = gower_dist(Y, centroids);  // n × K
  IntegerVector labels(n);
  for(int i = 0; i < n; ++i) {
    double best = Dk(i, 0);
    int   bidx = 0;
    for(int k = 1; k < K; ++k) {
      if (Dk(i, k) < best) {
        best = Dk(i, k);
        bidx = k;
      }
    }
    labels[i] = bidx;  // zero‐based state index
  }
  
  return labels;
}

static IntegerVector pam_cpp(const NumericMatrix& D, int K) {
  int n = D.nrow();
  IntegerVector medoids(K);
  std::vector<bool> is_med(n, false);
  // choose first medoid as row 0
  medoids[0] = 1;
  is_med[0] = true;
  std::vector<double> mindist(n);
  for (int i = 0; i < n; ++i) mindist[i] = D(i, 0);
  for (int k = 1; k < K; ++k) {
    int best_i = -1;
    double best_d = -1;
    for (int i = 0; i < n; ++i) {
      if (is_med[i]) continue;
      if (mindist[i] > best_d) {
        best_d = mindist[i];
        best_i = i;
      }
    }
    medoids[k] = best_i + 1;
    is_med[best_i] = true;
    for (int i = 0; i < n; ++i) {
      if (is_med[i]) continue;
      mindist[i] = std::min(mindist[i], D(i, best_i));
    }
  }
  return medoids;
}

// [[Rcpp::export]]
List robust_JM_COSA(NumericMatrix Y,
                    double zeta0,
                    double lambda,
                    int K,
                    Nullable<double> tol       = R_NilValue,
                    int n_init               = 10,
                    int n_outer              = 20,
                    double alpha            = 0.1,
                    bool verbose            = false,
                    int knn                 = 10,
                    double c                = 2.0,
                    Nullable<NumericMatrix> M = R_NilValue) {
  int TT = Y.nrow(), P = Y.ncol();
  bool use_tol = !tol.isNull();
  double tol_val = use_tol ? as<double>(tol) : 0.0;
  
  // Pre-allocate matrices/vectors
  NumericMatrix Gamma(K, K), W(K, P), W_old(K, P), DW, loss_by_state,
  Vdp(TT, K), Y_times(TT, P), Spk, wcd(K, P);
  IntegerVector s(TT), s_old;
  NumericVector v1, v2, v;
  double zeta, loss, loss_old, best_loss = R_PosInf;
  IntegerVector best_medoids;
  
  // Build transition penalty Γ
  for (int i = 0; i < K; ++i)
    for (int j = 0; j < K; ++j)
      Gamma(i, j) = (i != j ? lambda : 0.0);
  
  // Multi-init
  for (int init = 0; init < n_init; ++init) {
    W.fill(1.0 / P);
    W_old = clone(W);
    zeta = zeta0;
    loss_old = 1e10;
    
    s     = initialize_states(Y, K);
    s_old = clone(s);
    
    for (int outer = 0; outer < n_outer; ++outer) {
      // local scales
      for (int t = 0; t < TT; ++t) {
        int si = s[t] - 1;
        for (int p = 0; p < P; ++p)
          Y_times(t,p) = Y(t,p) * W(si,p);
      }
      if (M.isNotNull()) {
        v1 = Function("v_1")(Y_times, knn, c, M.get());
        v2 = Function("v_1")(Y,        knn, c, M.get());
      } else {
        v1 = Function("v_1")(Y_times, knn, c);
        v2 = Function("v_1")(Y,        knn, c);
      }
      v = pmin(v1, v2);
      
      // weighted distances + pam_cpp
      DW = weight_inv_exp_dist(Y_times, s, W, zeta);
      IntegerVector medoids = pam_cpp(DW, K);
      
      // loss by state
      loss_by_state = NumericMatrix(TT, K);
      for (int t = 0; t < TT; ++t)
        for (int k = 0; k < K; ++k)
          loss_by_state(t,k) = DW(t, medoids[k] - 1);
      
      // DP forward
      for (int k = 0; k < K; ++k)
        Vdp(TT-1,k) = loss_by_state(TT-1,k);
      for (int t = TT-2; t >= 0; --t) {
        for (int j = 0; j < K; ++j) {
          double best = R_PosInf;
          for (int i2 = 0; i2 < K; ++i2) {
            double cand = Vdp(t+1,i2) + Gamma(i2,j);
            if (cand < best) best = cand;
          }
          Vdp(t,j) = loss_by_state(t,j) + best;
        }
      }
      // backtrack
      {
        double bestv = Vdp(0,0); int bestk = 0;
        for (int k = 1; k < K; ++k)
          if (Vdp(0,k) < bestv) bestv = Vdp(0,k), bestk = k;
          s[0] = bestk + 1;
          loss  = bestv;
      }
      for (int t = 1; t < TT; ++t) {
        int prev = s[t-1] - 1;
        double bestv = R_PosInf; int bestk = 0;
        for (int k = 0; k < K; ++k) {
          double cand = Vdp(t,k) + Gamma(prev,k);
          if (cand < bestv) bestv = cand, bestk = k;
        }
        s[t] = bestk + 1;
      }
      
      // ensure K distinct
      if ((int)std::set<int>(s.begin(), s.end()).size() < K) {
        s = s_old;
        break;
      }
      s_old = clone(s);
      
      // loss convergence
      if (use_tol && (loss_old - loss < tol_val)) break;
      loss_old = loss;
      
      // update W
      Spk = WCD(s, Y_times, K);
      for (int ki = 0; ki < K; ++ki)
        for (int pi = 0; pi < P; ++pi)
          wcd(ki,pi) = std::exp(-Spk(ki,pi)/zeta0);
      for (int ki = 0; ki < K; ++ki) {
        double rs = 0;
        for (int pi = 0; pi < P; ++pi) rs += wcd(ki,pi);
        rs = (rs == 0 ? 1.0 : rs);
        for (int pi = 0; pi < P; ++pi)
          W(ki,pi) = wcd(ki,pi)/rs;
      }
      
      // W convergence
      double epsW = 0;
      for (int ki = 0; ki < K; ++ki)
        for (int pi = 0; pi < P; ++pi) {
          double d = W(ki,pi) - W_old(ki,pi);
          epsW += d*d;
        }
        epsW /= (K * P);
      if (use_tol && epsW < tol_val) break;
      W_old = clone(W);
      
      // bump zeta
      zeta += alpha * zeta0;
      
      if (verbose) {
        Rcout << "init=" << init+1
              << " outer=" << outer+1
              << " loss=" << loss
              << " epsW=" << epsW
              << " zeta=" << zeta << "\n";
      }
      
      if (loss < best_loss) {
        best_loss    = loss;
        best_medoids = medoids;
      }
    }
  }
  
  return List::create(
    _["W"]       = W,
    _["s"]       = s,
    _["medoids"] = best_medoids,
    _["v"]       = v,
    _["loss"]    = best_loss
  );
}

