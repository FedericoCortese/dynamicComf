// File: cpp_optimizer.cpp
#include <Rcpp.h>
#include <nlopt.hpp>
using namespace Rcpp;

// A little struct to carry data into the eval function
struct OptData {
  const std::vector<double> g;
  const std::vector<double> s_t1;
  double lambda, m;
};

// Objective+constraint definitions
double obj_fn(const std::vector<double>& x,
              std::vector<double>& grad,
              void* data_) {
  auto& d = *reinterpret_cast<OptData*>(data_);
  int K = x.size();
  double val=0, pen=0;
  for(int i=0;i<K;i++){
    val += std::pow(x[i], d.m) * d.g[i];
    double diff = d.s_t1[i] - x[i];
    pen += diff*diff;
    if (!grad.empty()) {
      // gradient = ∂/∂xᵢ [ sᵐ g + λ (sₜ₋₁ - s)² ]
      grad[i] = d.m * std::pow(x[i], d.m - 1)*d.g[i]
      - 2*d.lambda*diff;
    }
  }
  return val + d.lambda * pen;
}

// Equality constraint: sum(x) == 1
double eq_constr(const std::vector<double>& x,
                 std::vector<double>& grad,
                 void* /*data*/) {
  double s = std::accumulate(x.begin(), x.end(), 0.0);
  if (!grad.empty()) {
    for (size_t i = 0; i < x.size(); ++i)
      grad[i] = 1.0;
  }
  return s - 1.0;
}

// [[Rcpp::export]]
List optimize_cpp(NumericVector init_r,
                  NumericVector g_r,
                  NumericVector s_t1_r,
                  double lambda,
                  double m) {
  int K = init_r.size();
  // copy into std::vector
  std::vector<double> x(init_r.begin(), init_r.end());
  OptData data{ std::vector<double>(g_r.begin(), g_r.end()),
                std::vector<double>(s_t1_r.begin(), s_t1_r.end()),
                lambda, m };
  
  // Choose an algorithm (here MMA = Method of Moving Asymptotes)
  nlopt::opt opt(nlopt::LD_MMA, K);
  opt.set_min_objective(obj_fn, &data);
  opt.add_equality_constraint(eq_constr, nullptr, 1e-8);
  
  // bounds 0 ≤ xᵢ ≤ 1 (you can skip UB if you only need LB)
  std::vector<double> lb(K, 0.0), ub(K, 1.0);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  
  // stopping criteria
  opt.set_xtol_rel(1e-6);
  opt.set_maxeval(1e4);
  
  double minf;
  nlopt::result result = opt.optimize(x, minf);
  
  return List::create(
    _["par"]   = NumericVector(x.begin(), x.end()),
    _["value"] = minf,
    _["status"]= int(result)
  );
}
