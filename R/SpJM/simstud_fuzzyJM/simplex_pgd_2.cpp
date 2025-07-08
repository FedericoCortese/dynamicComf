// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Compute objective for 1T: sum(s^m * g) + lambda * sum(|s_t1 - s|^p_exp)
double eval_obj_1T(const NumericVector& s,
                   const NumericVector& g,
                   const NumericVector& s_t1,
                   double lambda,
                   double m,
                   int p_exp) {
  int K = s.size();
  double val = 0.0;
  double pen = 0.0;
  for (int i = 0; i < K; ++i) {
    val += std::pow(s[i], m) * g[i];
    double d = s_t1[i] - s[i];
    pen += std::pow(std::abs(d), p_exp);
  }
  return val + lambda * pen;
}

// Compute objective for 2T: sum(s^m * g) + lambda * (sum(|s_prec - s|^p_exp) + sum(|s_succ - s|^p_exp))
double eval_obj_2T(const NumericVector& s,
                   const NumericVector& g,
                   const NumericVector& s_prec,
                   const NumericVector& s_succ,
                   double lambda,
                   double m,
                   int p_exp) {
  int K = s.size();
  double val = 0.0;
  double pen = 0.0;
  for (int i = 0; i < K; ++i) {
    val += std::pow(s[i], m) * g[i];
    double d1 = s_prec[i] - s[i];
    double d2 = s_succ[i] - s[i];
    pen += std::pow(std::abs(d1), p_exp) + std::pow(std::abs(d2), p_exp);
  }
  return val + lambda * pen;
}

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Helper: project v onto the probability simplex { x >= 0, sum x = 1 }
NumericVector proj_simplex(const NumericVector& v) {
  int n = v.size();
  NumericVector u = clone(v);
  std::sort(u.begin(), u.end(), std::greater<double>());
  NumericVector css(n);
  css[0] = u[0];
  for (int i = 1; i < n; ++i) css[i] = css[i - 1] + u[i];
  double theta = 0;
  for (int i = n - 1; i >= 0; --i) {
    double t = (css[i] - 1.0) / (i + 1);
    if (u[i] > t) { theta = t; break; }
  }
  NumericVector w(n);
  for (int i = 0; i < n; ++i) w[i] = std::max(v[i] - theta, 0.0);
  return w;
}

// optimize single time-point with future neighbor s_t1
// Added time_norm: 1 = original squared-L2 penalty; 2 = L2-norm penalty
// [[Rcpp::export]]
List optimize_pgd_1T(NumericVector init,
                     NumericVector g,
                     NumericVector s_t1,
                     double lambda,
                     double m,
                     int time_norm = 1,
                     int max_iter = 1000,
                     double alpha    = 1e-2,
                     double tol      = 1e-8) {
  int K = init.size();
  NumericVector s = clone(init);
  for (int iter = 0; iter < max_iter; ++iter) {
    // Precompute norm_diff if needed
    double norm_diff = 0.0;
    if (time_norm == 2) {
      for (int i = 0; i < K; ++i) {
        double d = s_t1[i] - s[i]; norm_diff += d * d;
      }
      norm_diff = std::sqrt(norm_diff);
    }
    NumericVector grad(K);
    for (int i = 0; i < K; ++i) {
      // Data term gradient
      grad[i] = m * std::pow(s[i], m - 1) * g[i];
      // Temporal penalty gradient
      double d = s_t1[i] - s[i];
      if (time_norm == 1) {
        grad[i] += -2.0 * lambda * d;
      } else {
        if (norm_diff > 0) grad[i] += -lambda * (d / norm_diff);
      }
    }
    NumericVector s_new = proj_simplex(s - alpha * grad);
    if (max(abs(s_new - s)) < tol) { s = s_new; break; }
    s = s_new;
  }
  // Recompute objective components
  double val = 0.0, pen = 0.0, norm_diff = 0.0;
  for (int i = 0; i < K; ++i) {
    double d = s_t1[i] - s[i];
    val += std::pow(s[i], m) * g[i];
    if (time_norm == 1) {
      pen += d * d;
    } else {
      norm_diff += d * d;
    }
  }
  if (time_norm == 2) norm_diff = std::sqrt(norm_diff);
  double obj = val + (time_norm == 1 ? lambda * pen : lambda * norm_diff);
  return List::create(_["par"] = s, _["value"] = obj);
}

// optimize with both neighbors s_prec and s_succ
// time_norm: 1 = original squared-L2; 2 = L2-norm for neighbor differences
// [[Rcpp::export]]
List optimize_pgd_2T(NumericVector init,
                     NumericVector g,
                     NumericVector s_prec,
                     NumericVector s_succ,
                     double lambda,
                     double m,
                     int time_norm = 1,
                     int max_iter = 1000,
                     double alpha    = 1e-2,
                     double tol      = 1e-8) {
  int K = init.size();
  NumericVector s = clone(init);
  for (int iter = 0; iter < max_iter; ++iter) {
    double norm1 = 0.0, norm2 = 0.0;
    if (time_norm == 2) {
      for (int i = 0; i < K; ++i) {
        double d1 = s_prec[i] - s[i]; norm1 += d1 * d1;
        double d2 = s_succ[i] - s[i]; norm2 += d2 * d2;
      }
      norm1 = std::sqrt(norm1);
      norm2 = std::sqrt(norm2);
    }
    NumericVector grad(K);
    for (int i = 0; i < K; ++i) {
      grad[i] = m * std::pow(s[i], m - 1) * g[i];
      double d1 = s_prec[i] - s[i];
      if (time_norm == 1) {
        grad[i] += -2.0 * lambda * d1;
      } else if (norm1 > 0) {
        grad[i] += -lambda * (d1 / norm1);
      }
      double d2 = s_succ[i] - s[i];
      if (time_norm == 1) {
        grad[i] += -2.0 * lambda * d2;
      } else if (norm2 > 0) {
        grad[i] += -lambda * (d2 / norm2);
      }
    }
    NumericVector s_new = proj_simplex(s - alpha * grad);
    if (max(abs(s_new - s)) < tol) { s = s_new; break; }
    s = s_new;
  }
  // Recompute objective
  double val = 0.0, pen1 = 0.0, pen2 = 0.0;
  for (int i = 0; i < K; ++i) val += std::pow(s[i], m) * g[i];
  if (time_norm == 1) {
    for (int i = 0; i < K; ++i) {
      double d1 = s_prec[i] - s[i]; pen1 += d1 * d1;
      double d2 = s_succ[i] - s[i]; pen2 += d2 * d2;
    }
    val += lambda * (pen1 + pen2);
  } else {
    double norm1 = 0.0, norm2 = 0.0;
    for (int i = 0; i < K; ++i) {
      double d1 = s_prec[i] - s[i]; norm1 += d1 * d1;
      double d2 = s_succ[i] - s[i]; norm2 += d2 * d2;
    }
    val += lambda * (std::sqrt(norm1) + std::sqrt(norm2));
  }
  return List::create(_["par"] = s, _["value"] = val);
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

