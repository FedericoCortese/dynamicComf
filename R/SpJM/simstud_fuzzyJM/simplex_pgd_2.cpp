// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// ——— helper: project v onto simplex {x ≥ 0, ∑x = 1} ———
NumericVector proj_simplex(const NumericVector& v) {
  int n = v.size();
  NumericVector u = clone(v);
  std::sort(u.begin(), u.end(), std::greater<double>());
  NumericVector css(n);
  css[0] = u[0];
  for (int i = 1; i < n; ++i) css[i] = css[i - 1] + u[i];
  double rho = 0, theta = 0;
  for (int i = n - 1; i >= 0; --i) {
    double t = (css[i] - 1.0) / (i + 1);
    if (u[i] > t) { rho = i; theta = t; break; }
  }
  NumericVector w(n);
  for (int i = 0; i < n; ++i) w[i] = std::max(v[i] - theta, 0.0);
  return w;
}

// Compute objective for the 1T case: sum(s^m * g) + lambda * sum((s_t_1 - s)^2)
double eval_obj_1T(const NumericVector& s,
                   const NumericVector& g,
                   const NumericVector& s_t1,
                   double lambda,
                   double m) {
  int K = s.size();
  double val = 0.0, pen = 0.0;
  for (int i = 0; i < K; ++i) {
    val += std::pow(s[i], m) * g[i];
    double diff = s_t1[i] - s[i];
    pen += diff * diff;
  }
  return val + lambda * pen;
}

// Compute objective for the 2T case: sum(s^m * g) + lambda*(sum((s_t_prec - s)^2) + sum((s_t_succ - s)^2))
double eval_obj_2T(const NumericVector& s,
                   const NumericVector& g,
                   const NumericVector& s_prec,
                   const NumericVector& s_succ,
                   double lambda,
                   double m) {
  int K = s.size();
  double val = 0.0, pen = 0.0;
  for (int i = 0; i < K; ++i) {
    val += std::pow(s[i], m) * g[i];
    double d1 = s_prec[i] - s[i];
    double d2 = s_succ[i] - s[i];
    pen += d1 * d1 + d2 * d2;
  }
  return val + lambda * pen;
}

// [[Rcpp::export]]
List optimize_pgd_1T(NumericVector init,
                     NumericVector g,
                     NumericVector s_t1,
                     double lambda,
                     double m,
                     int max_iter = 1000,
                     double alpha    = 1e-2,
                     double tol      = 1e-8) {
  int K = init.size();
  NumericVector s = clone(init);
  for (int iter = 0; iter < max_iter; ++iter) {
    // gradient: m*s^(m-1)*g - 2*lambda*(s_t1 - s)
    NumericVector grad(K);
    for (int i = 0; i < K; ++i) {
      double diff = s_t1[i] - s[i];
      grad[i] = m * std::pow(s[i], m - 1) * g[i] - 2.0 * lambda * diff;
    }
    NumericVector s_new = proj_simplex(s - alpha * grad);
    if (max(abs(s_new - s)) < tol) {
      s = s_new;
      break;
    }
    s = s_new;
  }
  double obj = eval_obj_1T(s, g, s_t1, lambda, m);
  return List::create(
    _["par"]   = s,
    _["value"] = obj
  );
}

// [[Rcpp::export]]
List optimize_pgd_2T(NumericVector init,
                     NumericVector g,
                     NumericVector s_prec,
                     NumericVector s_succ,
                     double lambda,
                     double m,
                     int max_iter = 1000,
                     double alpha    = 1e-2,
                     double tol      = 1e-8) {
  int K = init.size();
  NumericVector s = clone(init);
  for (int iter = 0; iter < max_iter; ++iter) {
    // gradient: m*s^(m-1)*g - 2*lambda*((s_prec - s) + (s_succ - s))
    NumericVector grad(K);
    for (int i = 0; i < K; ++i) {
      double d1 = s_prec[i] - s[i];
      double d2 = s_succ[i] - s[i];
      grad[i] = m * std::pow(s[i], m - 1) * g[i] - 2.0 * lambda * (d1 + d2);
    }
    NumericVector s_new = proj_simplex(s - alpha * grad);
    if (max(abs(s_new - s)) < tol) {
      s = s_new;
      break;
    }
    s = s_new;
  }
  double obj = eval_obj_2T(s, g, s_prec, s_succ, lambda, m);
  return List::create(
    _["par"]   = s,
    _["value"] = obj
  );
}


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


// ——— fast weighted median ———
// [[Rcpp::export]]
double weighted_median(NumericVector x, NumericVector w) {
  int n = x.size();
  IntegerVector idx = seq(0, n - 1);
  std::sort(idx.begin(), idx.end(),
            [&](int i, int j){ return x[i] < x[j]; });
  double total = std::accumulate(w.begin(), w.end(), 0.0),
    half = total / 2, cum = 0;
  for (int r = 0; r < n; ++r) {
    cum += w[idx[r]];
    if (cum >= half) return x[idx[r]];
  }
  return x[idx[n-1]];
}

// ——— stub for initialize_states ———
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

// ——— inline loss (data + smoothness) ———
static inline double compute_loss(
    const NumericMatrix& S_pow,    // precomputed S^m
    const NumericMatrix& V,
    double lambda
) {
  int TT = V.nrow(), K = V.ncol();
  double L = 0.0;
  // data term
  for (int t = 0; t < TT; ++t)
    for (int k = 0; k < K; ++k)
      L += V(t, k) * S_pow(t, k);
  // smoothness
  for (int t = 0; t + 1 < TT; ++t)
    for (int k = 0; k < K; ++k) {
      double d = S_pow(t, k) - S_pow(t+1, k);
      L += lambda * d * d;
    }
    return L;
}

// [[Rcpp::export]]
List fit_sequence(const NumericMatrix& Y,
                       int K,
                       double lambda,
                       double m,
                       int n_init   = 10,
                       int max_iter = 1000,
                       double tol   = 1e-8,
                       bool verbose = false) {
  int TT = Y.nrow(), P = Y.ncol();
  
  // pre-allocate everything once
  NumericMatrix S       (TT, K),
  S_old   (TT, K),
  best_S  (TT, K),
  mu      (K, P),
  best_mu (K, P),
  V       (TT, K),
  S_pow   (TT, K);
  NumericVector init_vec(K), tmp_y(TT);
  double best_loss = R_PosInf;
  
  // uniform init vector
  for (int k = 0; k < K; ++k) init_vec[k] = 1.0 / K;
  
  // main n_init loop
  for (int run = 0; run < n_init; ++run) {
    // — 1) hard initialization via k-prototypes++ stub
    IntegerVector labels = initialize_states(Y, K);
    std::fill(S.begin(), S.end(), 0.0);
    for (int t = 0; t < TT; ++t)
      S(t, labels[t]) = 1.0;
    
    // — 2) initial μ via weighted medians
    //    need S_pow = S^m
    for (int t = 0; t < TT; ++t)
      for (int k = 0; k < K; ++k)
        S_pow(t, k) = std::pow(S(t, k), m);
    
    for (int k = 0; k < K; ++k) {
      // extract weights for this state
      NumericVector w = S_pow(_, k);
      for (int j = 0; j < P; ++j) {
        // build tmp_y = column j of Y
        for (int t = 0; t < TT; ++t) tmp_y[t] = Y(t, j);
        mu(k, j) = weighted_median(tmp_y, w);
      }
    }
    
    // — 3) initial distances & loss
    V = gower_dist(Y, mu);
    double loss_old = compute_loss(S_pow, V, lambda);
    S_old = clone(S);
    
    // — 4) PGD iterations
    for (int it = 0; it < max_iter; ++it) {
      // update S[0]
      {
        NumericVector par =
          optimize_pgd_1T(init_vec, V(0, _), S(1, _),
                          lambda, m)["par"];
        S(0, _) = par;
      }
      // update interior
      for (int t = 1; t + 1 < TT; ++t) {
        NumericVector par =
          optimize_pgd_2T(init_vec, V(t, _),
                          S(t-1, _), S(t+1, _),
                          lambda, m)["par"];
        S(t, _) = par;
      }
      // update last
      {
        NumericVector par =
          optimize_pgd_1T(init_vec, V(TT-1, _), S(TT-2, _),
                          lambda, m)["par"];
        S(TT-1, _) = par;
      }
      
      // recompute S^m
      for (int t = 0; t < TT; ++t)
        for (int k = 0; k < K; ++k)
          S_pow(t, k) = std::pow(S(t, k), m);
      
      // recompute μ
      for (int k = 0; k < K; ++k) {
        NumericVector w = S_pow(_, k);
        for (int j = 0; j < P; ++j) {
          for (int t = 0; t < TT; ++t) tmp_y[t] = Y(t, j);
          mu(k, j) = weighted_median(tmp_y, w);
        }
      }
      
      // recompute V & loss
      V = gower_dist(Y, mu);
      double loss = compute_loss(S_pow, V, lambda);
      if (verbose)
        Rcout << "init " << (run+1)
              << " iter " << (it+1)
              << " loss = " << loss << "\n";
        
        if ((loss_old - loss) < tol) break;
        loss_old = loss;
        S_old    = S;  // only shallow copy of values
    }
    
    // keep best
    if (loss_old < best_loss) {
      best_loss = loss_old;
      best_S    = clone(S);
      best_mu   = clone(mu);
    }
  }
  
  return List::create(
    _["S"]    = best_S,
    _["mu"]   = best_mu,
    _["loss"] = best_loss
  );
}


