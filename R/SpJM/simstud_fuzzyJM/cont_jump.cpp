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
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
using namespace Rcpp;

// Forward declarations
NumericMatrix discretize_prob_simplex(int K, double grid_size);
NumericMatrix gower_dist(const NumericMatrix& Y, const NumericMatrix& mu,
                         Nullable<IntegerVector> feat_type, std::string scale);

// helper: median of a vector
static double median_vec(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  int n = v.size();
  if(n == 0) return NA_REAL;
  return (n % 2 == 1) ? v[n/2] : 0.5 * (v[n/2 - 1] + v[n/2]);
}

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
               double grid_size = 0.05,
               bool verbose = false) {
  if(verbose) Rcout << "Starting cont_jump with " << n_init << " initializations\n";
  NumericMatrix Y_orig = clone(Y_in);
  int TT = Y_orig.nrow(), P = Y_orig.ncol();
  LogicalMatrix M_mask(TT, P);
  NumericVector colMean(P);
  for(int j=0;j<P;++j){ double sumj=0; int cntj=0;
  for(int t=0;t<TT;++t){ if(NumericVector::is_na(Y_orig(t,j))){ M_mask(t,j)=true; } else{ sumj+=Y_orig(t,j); cntj++; M_mask(t,j)=false;} }
  colMean[j] = cntj>0 ? sumj/cntj : 0.0;
  }
  for(int t=0;t<TT;++t) for(int j=0;j<P;++j) if(M_mask(t,j)) Y_orig(t,j)=colMean[j];
  
  NumericMatrix prob_vecs = discretize_prob_simplex(K, grid_size);
  int N = prob_vecs.nrow();
  NumericMatrix pairwise_l1(N,N);
  for(int i=0;i<N;++i) for(int j=0;j<N;++j){ double acc=0;
    for(int c=0;c<K;++c) acc+=std::abs(prob_vecs(i,c)-prob_vecs(j,c));
    pairwise_l1(i,j)=acc/2.0;
  }
  
  std::vector<List> all_runs; all_runs.reserve(n_init);
  
  for(int init_idx=0; init_idx<n_init; ++init_idx) {
    if(verbose) Rcout << " init "<< init_idx+1 <<"/"<<n_init<<"\n";
    NumericMatrix Y = clone(Y_orig);
    LogicalMatrix M = clone(M_mask);
    bool use_tol=false; double tol=0;
    if(tol_.isNotNull()){ tol=as<double>(tol_); use_tol=true; }
    IntegerVector s(TT);
    if(initial_states_.isNotNull()){
      IntegerVector init_s = initial_states_.get();
      if((int)init_s.size()!=TT) stop("initial_states must match rows");
      for(int t=0;t<TT;++t) s[t]=init_s[t]-1;
    } else {
      RNGScope scope;
      for(int t=0;t<TT;++t) s[t] = floor(R::runif(0, (double)K));
    }
    NumericMatrix mu(K,P);
    for(int k=0;k<K;++k){ int cnt=0; for(int j=0;j<P;++j) mu(k,j)=0;
    for(int t=0;t<TT;++t) if(s[t]==k){ cnt++; for(int j=0;j<P;++j) mu(k,j)+=Y(t,j);} 
    if(cnt>0) for(int j=0;j<P;++j) mu(k,j)/=cnt; else for(int j=0;j<P;++j) mu(k,j)=colMean[j];
    }
    
    NumericMatrix S(TT,K), S_old(TT,K), loss_mx(TT,N), values(TT,N);
    NumericVector assign(TT);
    double loss_old = R_PosInf, value_opt=R_PosInf;
    int iter_count=0;
    
    for(int iter=0; iter<max_iter; ++iter){ iter_count=iter+1;
      NumericMatrix loss_k(TT,K);
      for(int t=0;t<TT;++t) for(int k=0;k<K;++k){ double ss=0;
        for(int j=0;j<P;++j){ double d=Y(t,j)-mu(k,j); ss+=d*d;} loss_k(t,k)=0.5*sqrt(ss);
      }
      NumericMatrix jump_penalty_mx(N,N);
      for(int i=0;i<N;++i) for(int j=0;j<N;++j)
        jump_penalty_mx(i,j)=jump_penalty*pow(pairwise_l1(i,j),alpha);
      if(mode_loss){ NumericVector m_loss(N); 
        for(int i=0;i<N;++i){ double mx=-R_PosInf; for(int j=0;j<N;++j) mx = std::max(mx, -jump_penalty_mx(i,j));
        double se=0; for(int j=0;j<N;++j) se+=exp(-jump_penalty_mx(i,j)-mx);
        m_loss[i]=mx+log(se);
        }
        double off=m_loss[0]; for(int i=0;i<N;++i) m_loss[i]-=off;
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) jump_penalty_mx(i,j)+=m_loss[i];
      }
      for(int t=0;t<TT;++t) for(int i=0;i<N;++i){ double ac=0; for(int k=0;k<K;++k) ac+=loss_k(t,k)*prob_vecs(i,k); loss_mx(t,i)=ac; }
      for(int t=0;t<TT;++t) for(int i=0;i<N;++i) if(std::isnan(loss_mx(t,i))) loss_mx(t,i)=R_PosInf;
      for(int i=0;i<N;++i) values(0,i)=loss_mx(0,i);
      for(int t=1;t<TT;++t) for(int i=0;i<N;++i){ double bv=R_PosInf; for(int j=0;j<N;++j) bv=std::min(bv, values(t-1,j)+jump_penalty_mx(j,i)); values(t,i)=loss_mx(t,i)+bv; }
      double best=R_PosInf; for(int i=0;i<N;++i) if(values(TT-1,i)<best){ best=values(TT-1,i); assign[TT-1]=i; }
      value_opt=best;
      for(int k=0;k<K;++k) S(TT-1,k)=prob_vecs(assign[TT-1],k);
      for(int t=TT-2;t>=0;--t){ int nxt=assign[t+1]; double bv=R_PosInf; int bi=0;
      for(int i=0;i<N;++i){ double cnd=values(t,i)+jump_penalty_mx(i,nxt); if(cnd<bv){bv=cnd;bi=i;}} assign[t]=bi;
      for(int k=0;k<K;++k) S(t,k)=prob_vecs(bi,k);
      }
      NumericVector sum_prob(K); for(int k=0;k<K;++k){ sum_prob[k]=0; for(int j=0;j<P;++j) mu(k,j)=0; }
      for(int t=0;t<TT;++t) for(int k=0;k<K;++k){ sum_prob[k]+=S(t,k); for(int j=0;j<P;++j) mu(k,j)+=S(t,k)*Y(t,j);}      
      for(int k=0;k<K;++k) if(sum_prob[k]>0) for(int j=0;j<P;++j) mu(k,j)/=sum_prob[k];
      for(int t=0;t<TT;++t){ int hk=0; double mp=S(t,0); for(int k=1;k<K;++k) if(S(t,k)>mp){mp=S(t,k);hk=k;} for(int j=0;j<P;++j) if(M_mask(t,j)) Y(t,j)=mu(hk,j); }
      bool converged=true;
      if(use_tol){ double eps=loss_old-value_opt; if(eps<tol) converged=true; else converged=false; }
      else { converged=true; for(int t=0;t<TT && converged;++t) for(int k=0;k<K;++k) if(S(t,k)!=S_old(t,k)){ converged=false; break;} }
      if(converged) { if(verbose) Rcout<<"  converged at iter="<<iter_count<<" loss="<<value_opt<<"\n"; break; }
      loss_old=value_opt;
      for(int t=0;t<TT;++t) for(int k=0;k<K;++k) S_old(t,k)=S(t,k);
    }
    if(verbose) Rcout<<" end init "<<init_idx+1<<" loss="<<value_opt<<" iters="<<iter_count<<"\n";
    all_runs.push_back(List::create(Named("S")=S, Named("value_opt")=value_opt, Named("mu")=mu));
  }
  double best_loss=R_PosInf; int best_idx=0;
  for(int i=0;i<n_init;++i){ double v=all_runs[i]["value_opt"]; if(v<best_loss){best_loss=v;best_idx=i;} }
  if(verbose) Rcout<<"Best init="<<best_idx+1<<" loss="<<best_loss<<"\n";
  List best=all_runs[best_idx];
  return List::create(Named("best_S")=best["S"], Named("best_loss")=best["value_opt"], Named("best_mu")=best["mu"]);
}

