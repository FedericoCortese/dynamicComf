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
  for (int i = 1; i < n; ++i) {
    css[i] = css[i - 1] + u[i];
  }
  double rho = 0, theta = 0;
  for (int i = n - 1; i >= 0; --i) {
    double t = (css[i] - 1.0) / (i + 1);
    if (u[i] > t) {
      rho = i;
      theta = t;
      break;
    }
  }
  NumericVector w(n);
  for (int i = 0; i < n; ++i) {
    w[i] = std::max(v[i] - theta, 0.0);
  }
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