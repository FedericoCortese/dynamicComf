// cont_jump.cpp
// Compile with: Rcpp::sourceCpp("cont_jump.cpp")

#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------------------------------------
// Helper: generate all integer combinations of length n_c that sum to N
//------------------------------------------------------------------------------
static void generate_combinations(int n_c, int N, 
                                  std::vector< std::vector<int> > &out,
                                  std::vector<int> &current,
                                  int idx, int remaining) 
{
  if (idx == n_c - 1) {
    current[idx] = remaining;
    out.push_back(current);
    return;
  }
  for (int i = 0; i <= remaining; ++i) {
    current[idx] = i;
    generate_combinations(n_c, N, out, current, idx + 1, remaining - i);
  }
}

//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix discretize_prob_simplex(int n_c, double grid_size) {
  // Sample a grid on the n_c-dimensional probability simplex.
  int N = (int) std::round(1.0 / grid_size);
  if (N <= 0) {
    stop("grid_size must be a positive fraction <= 1");
  }
  
  // Generate all integer tuples of length n_c that sum to N
  std::vector< std::vector<int> > int_tuples;
  std::vector<int> current(n_c, 0);
  generate_combinations(n_c, N, int_tuples, current, 0, N);
  
  int M = (int) int_tuples.size();
  NumericMatrix simplex(M, n_c);
  
  // Fill in reverse order to match R code’s reverse indexing
  for (int i = 0; i < M; ++i) {
    int rev_i = M - 1 - i;
    for (int j = 0; j < n_c; ++j) {
      simplex(i, j) = double(int_tuples[rev_i][j]) / double(N);
    }
  }
  
  return simplex;
}

//------------------------------------------------------------------------------
// [[Rcpp::export]]
List cont_jump(const NumericMatrix &Y_in,
               int K,
               double jump_penalty = 1e-5,
               double alpha = 2.0,
               Nullable<IntegerVector> initial_states_ = R_NilValue,
               int max_iter = 10,
               int n_init = 10,
               Nullable<double> tol_ = R_NilValue,
               bool mode_loss = true,
               double grid_size = 0.05) 
{
  // Fit continuous jump model with 'n_init' random initializations.
  // Returns the best S, best_loss, and best_mu.
  
  // 1. Clone input and build missing mask M
  NumericMatrix Y_orig = clone(Y_in);
  int TT = Y_orig.nrow();
  int P  = Y_orig.ncol();
  LogicalMatrix M_mask(TT, P);
  
  // Compute column means (ignoring NA) for initial imputation
  NumericVector colMean(P, 0.0);
  for (int j = 0; j < P; ++j) {
    double sumj = 0.0;
    int cntj = 0;
    for (int t = 0; t < TT; ++t) {
      if (NumericVector::is_na(Y_orig(t, j))) {
        M_mask(t, j) = true;
      } else {
        sumj += Y_orig(t, j);
        cntj++;
        M_mask(t, j) = false;
      }
    }
    colMean[j] = (cntj > 0 ? sumj / double(cntj) : 0.0);
  }
  // Impute missing with column means
  for (int t = 0; t < TT; ++t) {
    for (int j = 0; j < P; ++j) {
      if (M_mask(t, j)) {
        Y_orig(t, j) = colMean[j];
      }
    }
  }
  
  // 2. Precompute simplex grid and pairwise distances
  NumericMatrix prob_vecs = discretize_prob_simplex(K, grid_size);
  int N = prob_vecs.nrow();
  NumericMatrix pairwise_l1(N, N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      double acc = 0.0;
      for (int c = 0; c < K; ++c) {
        acc += std::abs(prob_vecs(i, c) - prob_vecs(j, c));
      }
      pairwise_l1(i, j) = acc / 2.0;
    }
  }
  
  // 3. Container to store results of each initialization
  std::vector< List > all_runs;
  all_runs.reserve(n_init);
  
  // 4. Run n_init random initializations
  for (int init_idx = 0; init_idx < n_init; ++init_idx) {
    // 4a. Clone Y for this run
    NumericMatrix Y = clone(Y_orig);
    LogicalMatrix M = clone(M_mask);
    
    // 4b. Handle tol
    double tol = 0.0;
    bool use_tol = false;
    if (tol_.isNotNull()) {
      tol = as<double>(tol_);
      use_tol = true;
    }
    
    // 4c. Initialize state assignments s[t] in [0..K-1]
    IntegerVector s(TT);
    if (initial_states_.isNotNull()) {
      IntegerVector init_s = as<IntegerVector>(initial_states_);
      if ((int) init_s.size() != TT) {
        stop("initial_states must have length == nrow(Y)");
      }
      for (int t = 0; t < TT; ++t) {
        int st = init_s[t];
        if (st < 1 || st > K) stop("initial_states values must be in 1..K");
        s[t] = st - 1; // 0-based
      }
    } else {
      RNGScope scope;
      for (int t = 0; t < TT; ++t) {
        int rv = (int) std::floor(R::runif(0.0, (double)K));
        if (rv < 0) rv = 0;
        if (rv >= K) rv = K - 1;
        s[t] = rv;
      }
    }
    
    // 4d. Initialize mu[K x P]
    NumericMatrix mu(K, P);
    {
      for (int k = 0; k < K; ++k) {
        int count_k = 0;
        for (int j = 0; j < P; ++j) {
          mu(k, j) = 0.0;
        }
        for (int t = 0; t < TT; ++t) {
          if (s[t] == k) {
            count_k++;
            for (int j = 0; j < P; ++j) {
              mu(k, j) += Y(t, j);
            }
          }
        }
        if (count_k > 0) {
          for (int j = 0; j < P; ++j) {
            mu(k, j) /= double(count_k);
          }
        } else {
          // If no points in cluster, set mu[k,] to column means
          for (int j = 0; j < P; ++j) {
            mu(k, j) = colMean[j];
          }
        }
      }
    }
    
    // 4e. Prepare for EM‐DP iterations
    NumericMatrix S(TT, K), S_old(TT, K);
    double loss_old = R_PosInf;
    double value_opt = R_PosInf;
    
    NumericMatrix loss_mx(TT, N), values(TT, N);
    IntegerVector assign(TT);
    
    // 4f. EM iterations
    for (int iter = 0; iter < max_iter; ++iter) {
      // E-step: compute loss_k[t,k] = 0.5 * || Y[t,] - mu[k,] ||_2
      NumericMatrix loss_k(TT, K);
      for (int t = 0; t < TT; ++t) {
        for (int k = 0; k < K; ++k) {
          double sumsq = 0.0;
          for (int j = 0; j < P; ++j) {
            double diff = Y(t, j) - mu(k, j);
            sumsq += diff * diff;
          }
          loss_k(t, k) = 0.5 * std::sqrt(sumsq);
        }
      }
      // Build jump_penalty_mx = jump_penalty * (pairwise_l1 ^ alpha)
      NumericMatrix jump_penalty_mx(N, N);
      for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
          jump_penalty_mx(i, j) = jump_penalty * std::pow(pairwise_l1(i, j), alpha);
        }
      }
      // mode_loss adjustment
      if (mode_loss) {
        NumericVector m_loss(N, 0.0);
        for (int i = 0; i < N; ++i) {
          double row_max = R_NegInf;
          for (int j = 0; j < N; ++j) {
            double v = -jump_penalty_mx(i, j);
            if (v > row_max) row_max = v;
          }
          double sumexp = 0.0;
          for (int j = 0; j < N; ++j) {
            sumexp += std::exp(-jump_penalty_mx(i, j) - row_max);
          }
          m_loss[i] = row_max + std::log(sumexp);
        }
        double offset = m_loss[0];
        for (int i = 0; i < N; ++i) {
          m_loss[i] -= offset;
        }
        for (int i = 0; i < N; ++i) {
          for (int j = 0; j < N; ++j) {
            jump_penalty_mx(i, j) += m_loss[i];
          }
        }
      }
      // Compute loss_mx = loss_k %*% t(prob_vecs)
      for (int t = 0; t < TT; ++t) {
        for (int i = 0; i < N; ++i) {
          double acc = 0.0;
          for (int k = 0; k < K; ++k) {
            acc += loss_k(t, k) * prob_vecs(i, k);
          }
          loss_mx(t, i) = acc;
        }
      }
      // Replace NaN with Inf
      for (int t = 0; t < TT; ++t) {
        for (int i = 0; i < N; ++i) {
          if (std::isnan(loss_mx(t, i))) {
            loss_mx(t, i) = R_PosInf;
          }
        }
      }
      // DP: values[0, ] = loss_mx[0, ]
      for (int i = 0; i < N; ++i) {
        values(0, i) = loss_mx(0, i);
      }
      // Recursion
      for (int t = 1; t < TT; ++t) {
        for (int i = 0; i < N; ++i) {
          double bestv = R_PosInf;
          for (int j = 0; j < N; ++j) {
            double cand = values(t - 1, j) + jump_penalty_mx(j, i);
            if (cand < bestv) bestv = cand;
          }
          values(t, i) = loss_mx(t, i) + bestv;
        }
      }
      // Backtrack end
      double best_end = R_PosInf;
      for (int i = 0; i < N; ++i) {
        if (values(TT - 1, i) < best_end) {
          best_end    = values(TT - 1, i);
          assign[TT - 1] = i;
        }
      }
      value_opt = best_end;
      for (int k = 0; k < K; ++k) {
        S(TT - 1, k) = prob_vecs(assign[TT - 1], k);
      }
      // Traceback
      for (int t = TT - 2; t >= 0; --t) {
        int next_idx = assign[t + 1];
        double bestv = R_PosInf;
        int best_i   = 0;
        for (int i = 0; i < N; ++i) {
          double cand = values(t, i) + jump_penalty_mx(i, next_idx);
          if (cand < bestv) {
            bestv = cand;
            best_i = i;
          }
        }
        assign[t] = best_i;
        for (int k = 0; k < K; ++k) {
          S(t, k) = prob_vecs(best_i, k);
        }
      }
      // M-step: update mu
      NumericVector sum_prob(K, 0.0);
      for (int k = 0; k < K; ++k) {
        for (int j = 0; j < P; ++j) {
          mu(k, j) = 0.0;
        }
      }
      for (int t = 0; t < TT; ++t) {
        for (int k = 0; k < K; ++k) {
          sum_prob[k] += S(t, k);
          for (int j = 0; j < P; ++j) {
            mu(k, j) += S(t, k) * Y(t, j);
          }
        }
      }
      for (int k = 0; k < K; ++k) {
        if (sum_prob[k] > 0) {
          for (int j = 0; j < P; ++j) {
            mu(k, j) /= sum_prob[k];
          }
        }
      }
      // Re-impute missing in Y
      for (int t = 0; t < TT; ++t) {
        int hard_k = 0;
        double maxp = S(t, 0);
        for (int k = 1; k < K; ++k) {
          if (S(t, k) > maxp) {
            maxp = S(t, k);
            hard_k = k;
          }
        }
        for (int j = 0; j < P; ++j) {
          if (M(t, j)) {
            Y(t, j) = mu(hard_k, j);
          }
        }
      }
      // Check convergence
      if (use_tol) {
        double eps = loss_old - value_opt;
        if (eps < tol) break;
      } else {
        bool same = true;
        for (int t = 0; t < TT && same; ++t) {
          for (int k = 0; k < K; ++k) {
            if (S(t, k) != S_old(t, k)) {
              same = false;
              break;
            }
          }
        }
        if (same) break;
      }
      loss_old = value_opt;
      // Copy S to S_old
      for (int t = 0; t < TT; ++t) {
        for (int k = 0; k < K; ++k) {
          S_old(t, k) = S(t, k);
        }
      }
    } // end EM loop
    
    // Store this run's results
    all_runs.push_back(List::create(
        Named("S")         = S,
        Named("value_opt") = value_opt,
        Named("mu")        = mu
    ));
  } // end n_init loop
  
  // 5. Select best run (minimum value_opt)
  double best_loss = R_PosInf;
  int best_idx = 0;
  for (int i = 0; i < n_init; ++i) {
    double val = as<double>(all_runs[i]["value_opt"]);
    if (val < best_loss) {
      best_loss = val;
      best_idx = i;
    }
  }
  List best_run = all_runs[best_idx];
  
  return List::create(
    Named("best_S")    = best_run["S"],
                                 Named("best_loss") = best_run["value_opt"],
                                                              Named("best_mu")   = best_run["mu"]
  );
}
