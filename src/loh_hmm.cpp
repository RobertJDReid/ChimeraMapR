// =============================================================================
//  loh_hmm.cpp
//
//  Rcpp port of the beta-binomial EM + Viterbi HMM used by compute_loh_map()
//  in chimera_functions.R. Same model, same numerics — moved out of R only
//  because the EM's per-iteration optimize() calls and the per-position
//  Viterbi for-loop scale with SNP count and were the hot path for large
//  genomes. Compiled on demand via Rcpp::sourceCpp() from chimera_functions.R.
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <array>
#include <vector>

using namespace Rcpp;

static inline double lbeta_(double a, double b) {
  return R::lgammafn(a) + R::lgammafn(b) - R::lgammafn(a + b);
}

static inline double lchoose_(double n, double k) {
  return R::lgammafn(n + 1.0) - R::lgammafn(k + 1.0) - R::lgammafn(n - k + 1.0);
}

static inline double dbetabinom1(double k, double n, double mu, double theta,
                                  double lbab) {
  double alpha = mu * theta;
  double beta_ = (1.0 - mu) * theta;
  double val = std::exp(lchoose_(n, k) + lbeta_(k + alpha, n - k + beta_) - lbab);
  return std::isfinite(val) ? val : 0.0;
}

// [[Rcpp::export]]
NumericVector dbetabinom_cpp(NumericVector k, NumericVector N, double mu, double theta) {
  int n = k.size();
  NumericVector out(n);
  double alpha = mu * theta;
  double beta_ = (1.0 - mu) * theta;
  double lbab  = lbeta_(alpha, beta_);
  for (int i = 0; i < n; i++) out[i] = dbetabinom1(k[i], N[i], mu, theta, lbab);
  return out;
}

// Weighted log-likelihood of the HET beta-binomial component (mu = 0.5),
// as a function of the overdispersion parameter theta. Mirrors wll_theta()
// in the R implementation.
static double wll_theta(double theta, const NumericVector &k, const NumericVector &N,
                        const std::vector<double> &w) {
  double alpha = 0.5 * theta, beta_ = 0.5 * theta;
  double lbab  = lbeta_(alpha, beta_);
  double total = 0.0;
  int n = k.size();
  for (int i = 0; i < n; i++) {
    double val = dbetabinom1(k[i], N[i], 0.5, theta, lbab);
    if (val <= 0) val = 1e-300;
    total += w[i] * std::log(val);
  }
  return total;
}

// Golden-section search maximizing wll_theta() over [lo, hi]. Replaces R's
// optimize() — same idea (no derivatives needed), just without the R-call
// overhead of re-entering this from every EM iteration.
static double maximize_theta(const NumericVector &k, const NumericVector &N,
                              const std::vector<double> &w,
                              double lo, double hi, double tol) {
  const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
  double a = lo, b = hi;
  double c = b - gr * (b - a);
  double d = a + gr * (b - a);
  double fc = wll_theta(c, k, N, w);
  double fd = wll_theta(d, k, N, w);
  while (std::fabs(b - a) > tol) {
    if (fc > fd) {
      b = d; d = c; fd = fc;
      c = b - gr * (b - a);
      fc = wll_theta(c, k, N, w);
    } else {
      a = c; c = d; fc = fd;
      d = a + gr * (b - a);
      fd = wll_theta(d, k, N, w);
    }
  }
  return (a + b) / 2.0;
}

// 3-state mixture: HOM_REF (binomial, error rate eps), DIP_HET_0.5
// (beta-binomial, mu=0.5), HOM_ALT (binomial, error rate 1-eps).
// Mirrors fit_ab_mixture() in compute_loh_map().
// [[Rcpp::export]]
List fit_ab_mixture_cpp(NumericVector k, NumericVector N,
                        double theta_init = 80.0, double eps_init = 0.02,
                        int max_iter = 200, double tol = 1e-5) {
  int n = k.size();
  double pi0 = 1.0 / 3.0, pi1 = 1.0 / 3.0, pi2 = 1.0 / 3.0;
  double eps   = eps_init;
  double theta = theta_init;

  std::vector<double> r0(n), r1(n), r2(n);

  for (int iter = 0; iter < max_iter; iter++) {
    double old_pi0 = pi0, old_pi1 = pi1, old_pi2 = pi2, old_eps = eps, old_theta = theta;

    double eps_lo = std::min(std::max(eps, 1e-9), 1.0 - 1e-9);
    double eps_hi = std::min(std::max(1.0 - eps, 1e-9), 1.0 - 1e-9);
    double alpha = 0.5 * theta, beta_ = 0.5 * theta;
    double lbab  = lbeta_(alpha, beta_);

    for (int i = 0; i < n; i++) {
      double d0 = R::dbinom(k[i], N[i], eps_lo, 0);
      double d1 = dbetabinom1(k[i], N[i], 0.5, theta, lbab);
      double d2 = R::dbinom(k[i], N[i], eps_hi, 0);

      double w0 = pi0 * d0, w1 = pi1 * d1, w2 = pi2 * d2;
      double rs = w0 + w1 + w2;
      if (rs <= 0) rs = 1e-300;
      r0[i] = w0 / rs; r1[i] = w1 / rs; r2[i] = w2 / rs;
    }

    double sum0 = 0, sum1 = 0, sum2 = 0;
    for (int i = 0; i < n; i++) { sum0 += r0[i]; sum1 += r1[i]; sum2 += r2[i]; }
    pi0 = sum0 / n; pi1 = sum1 / n; pi2 = sum2 / n;

    double num = 0, den = 0;
    for (int i = 0; i < n; i++) {
      num += r0[i] * (k[i] / N[i]) + r2[i] * (1.0 - k[i] / N[i]);
      den += r0[i] + r2[i];
    }
    if (den > 0) eps = std::min(std::max(num / den, 1e-4), 0.15);

    if (sum1 >= 1.0) theta = maximize_theta(k, N, r1, 5.0, 5000.0, 1e-4);

    double dpi0   = std::fabs(pi0 - old_pi0)     / (std::fabs(old_pi0) + 1.0);
    double dpi1   = std::fabs(pi1 - old_pi1)     / (std::fabs(old_pi1) + 1.0);
    double dpi2   = std::fabs(pi2 - old_pi2)     / (std::fabs(old_pi2) + 1.0);
    double deps   = std::fabs(eps - old_eps)     / (std::fabs(old_eps) + 1.0);
    double dtheta = std::fabs(theta - old_theta) / (std::fabs(old_theta) + 1.0);
    double maxdiff = std::max({dpi0, dpi1, dpi2, deps, dtheta});
    if (maxdiff < tol) break;
  }

  return List::create(
    Named("mixing_proportions") = NumericVector::create(pi0, pi1, pi2),
    Named("error_rate")         = eps,
    Named("theta")              = theta
  );
}

// Per-chromosome Viterbi decoding over the 3-state mixture fit above.
// `chrom` must be pre-sorted (chrom, then pos) — chromosome runs are taken
// as contiguous blocks, matching the setorder() done in compute_loh_map()
// before this is called.
// [[Rcpp::export]]
CharacterVector viterbi_segment_cpp(NumericVector k, NumericVector N, CharacterVector chrom,
                                    NumericVector mixing_proportions,
                                    double error_rate, double theta, double trans_stay) {
  int n = k.size();
  const int S = 3;

  double trans_off = (1.0 - trans_stay) / (S - 1);
  double A_log[S][S];
  for (int i = 0; i < S; i++)
    for (int j = 0; j < S; j++)
      A_log[i][j] = (i == j) ? std::log(trans_stay) : std::log(trans_off);

  double eps_lo = std::min(std::max(error_rate, 1e-9), 1.0 - 1e-9);
  double eps_hi = std::min(std::max(1.0 - error_rate, 1e-9), 1.0 - 1e-9);
  double alpha = 0.5 * theta, beta_ = 0.5 * theta;
  double lbab  = lbeta_(alpha, beta_);

  double log_pi[S];
  for (int s = 0; s < S; s++) log_pi[s] = std::log(mixing_proportions[s]);

  std::vector<std::array<double, S>> log_emit(n);
  for (int i = 0; i < n; i++) {
    double e0 = R::dbinom(k[i], N[i], eps_lo, 0);
    double e1 = dbetabinom1(k[i], N[i], 0.5, theta, lbab);
    double e2 = R::dbinom(k[i], N[i], eps_hi, 0);
    log_emit[i][0] = std::log(std::max(e0, 1e-300));
    log_emit[i][1] = std::log(std::max(e1, 1e-300));
    log_emit[i][2] = std::log(std::max(e2, 1e-300));
  }

  CharacterVector labels = CharacterVector::create("HOM_REF", "DIP_HET_0.5", "HOM_ALT");
  CharacterVector assignments(n);

  int i = 0;
  while (i < n) {
    std::string cur = as<std::string>(chrom[i]);
    int j = i + 1;
    while (j < n && as<std::string>(chrom[j]) == cur) j++;
    int T_ = j - i;

    std::vector<std::array<double, S>> delta(T_);
    std::vector<std::array<int, S>>    psi(T_);

    for (int s = 0; s < S; s++) delta[0][s] = log_pi[s] + log_emit[i][s];

    for (int t = 1; t < T_; t++) {
      for (int s2 = 0; s2 < S; s2++) {
        double best = R_NegInf;
        int bestArg = 0;
        for (int s1 = 0; s1 < S; s1++) {
          double score = delta[t - 1][s1] + A_log[s1][s2];
          if (score > best) { best = score; bestArg = s1; }
        }
        psi[t][s2]   = bestArg;
        delta[t][s2] = best + log_emit[i + t][s2];
      }
    }

    int path_T = 0;
    double best = delta[T_ - 1][0];
    for (int s = 1; s < S; s++)
      if (delta[T_ - 1][s] > best) { best = delta[T_ - 1][s]; path_T = s; }

    std::vector<int> path(T_);
    path[T_ - 1] = path_T;
    for (int t = T_ - 2; t >= 0; t--) path[t] = psi[t + 1][path[t + 1]];

    for (int t = 0; t < T_; t++) assignments[i + t] = labels[path[t]];

    i = j;
  }

  return assignments;
}
