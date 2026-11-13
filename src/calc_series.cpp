#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

//' Evaluates Moschopoulos Truncated Series CDF / Survival Function
//'
//' Fast C++ vectorized evaluation of the CDF or survival function of a sum of Gammas
//' using the precomputed Moschopoulos mixture weights with Poisson-bracketed active windows.
//'
//' @param x_in NumericVector of evaluation points.
//' @param sumalpha Numeric, sum of shapes.
//' @param tau0 Numeric, base rate (maximum rate of the sorted components).
//' @param i0 Integer, 1-based start index of the trimmed mixture weights.
//' @param w NumericVector of mixture weights.
//' @param cw_lower NumericVector, prefix cumulative sums of w.
//' @param cw_upper NumericVector, suffix cumulative sums of w.
//' @param lower_tail Logical, if TRUE, returns P(X <= x); if FALSE, returns P(X > x).
//' @param log_p Logical, if TRUE, returns log probabilities.
//' @return NumericVector of evaluated probabilities.
//' @export
// [[Rcpp::export]]
NumericVector psumgamma_series_cpp(NumericVector x_in, double sumalpha, double tau0, int i0,
                                  NumericVector w, NumericVector cw_lower, NumericVector cw_upper,
                                  bool lower_tail = true, bool log_p = false) {
  const int nx = x_in.size();
  const int kw = w.size();
  NumericVector out(nx);
  const double scale0 = 1.0 / tau0;
  const double* p_w = REAL(w);
  const double* p_cw_low = REAL(cw_lower);
  const double* p_cw_up = REAL(cw_upper);
  
  for (int idx = 0; idx < nx; ++idx) {
    double xi = x_in[idx];
    if (xi <= 0.0) {
      double val = lower_tail ? 0.0 : 1.0;
      out[idx] = log_p ? (val > 0.0 ? 0.0 : -R_PosInf) : val;
      continue;
    }
    if (std::isinf(xi)) {
      double val = lower_tail ? 1.0 : 0.0;
      out[idx] = log_p ? (val > 0.0 ? 0.0 : -R_PosInf) : val;
      continue;
    }
    
    double lambda = xi * tau0;
    double a_min = R::qpois(1e-16, lambda, 1, 0);
    double a_max = R::qpois(1e-16, lambda, 0, 0);
    
    int jm = std::min(kw - 1, std::max(0, static_cast<int>(std::floor(a_min - sumalpha)) - 4 - (i0 - 1)));
    int jM = std::max(jm, std::min(kw - 1, static_cast<int>(std::ceil(a_max - sumalpha)) + 6 - (i0 - 1)));
    
    double p_base = lower_tail ? p_cw_low[jm] : p_cw_up[jM + 1];
    double p_sum = 0.0;
    for (int j = jm; j <= jM; ++j) {
      double ak = sumalpha + (i0 - 1) + j;
      p_sum += p_w[j] * R::pgamma(xi, ak, scale0, lower_tail ? 1 : 0, 0);
    }
    double total_p = p_base + p_sum;
    if (log_p) {
      out[idx] = (total_p > 0.0) ? std::log(total_p) : -R_PosInf;
    } else {
      out[idx] = total_p;
    }
  }
  return out;
}

//' Evaluates Moschopoulos Truncated Series Probability Density Function
//'
//' Fast C++ vectorized evaluation of the density function of a sum of Gammas
//' using the precomputed Moschopoulos mixture weights with Poisson-bracketed active windows.
//'
//' @param x_in NumericVector of evaluation points.
//' @param sumalpha Numeric, sum of shapes.
//' @param tau0 Numeric, base rate (maximum rate of the sorted components).
//' @param i0 Integer, 1-based start index of the trimmed mixture weights.
//' @param w NumericVector of mixture weights.
//' @param log_d Logical, if TRUE, returns log densities.
//' @return NumericVector of evaluated densities.
//' @export
// [[Rcpp::export]]
NumericVector dsumgamma_series_cpp(NumericVector x_in, double sumalpha, double tau0, int i0,
                                  NumericVector w, bool log_d = false) {
  const int nx = x_in.size();
  const int kw = w.size();
  NumericVector out(nx);
  const double scale0 = 1.0 / tau0;
  const double* p_w = REAL(w);
  
  for (int idx = 0; idx < nx; ++idx) {
    double xi = x_in[idx];
    if (xi <= 0.0) {
      out[idx] = log_d ? -R_PosInf : 0.0;
      continue;
    }
    if (std::isinf(xi)) {
      out[idx] = log_d ? -R_PosInf : 0.0;
      continue;
    }
    
    double lambda = xi * tau0;
    double a_min = R::qpois(1e-16, lambda, 1, 0);
    double a_max = R::qpois(1e-16, lambda, 0, 0);
    
    int jm = std::min(kw - 1, std::max(0, static_cast<int>(std::floor(a_min - sumalpha)) - 4 - (i0 - 1)));
    int jM = std::max(jm, std::min(kw - 1, static_cast<int>(std::ceil(a_max - sumalpha)) + 6 - (i0 - 1)));
    
    double d_sum = 0.0;
    for (int j = jm; j <= jM; ++j) {
      double ak = sumalpha + (i0 - 1) + j;
      d_sum += p_w[j] * R::dgamma(xi, ak, scale0, 0);
    }
    
    if (log_d) {
      out[idx] = (d_sum > 0.0) ? std::log(d_sum) : -R_PosInf;
    } else {
      out[idx] = d_sum;
    }
  }
  return out;
}

//' Exact Random Variate Generation from Moschopoulos Mixture
//'
//' Draws n random variates from the sum of Gammas via exact two-stage categorical
//' mixture sampling.
//'
//' @param n Integer, number of variates to sample.
//' @param sumalpha Numeric, sum of shapes.
//' @param tau0 Numeric, base rate.
//' @param i0 Integer, 1-based start index of weights.
//' @param cw_lower NumericVector, prefix cumulative sums of weights (length kw + 1).
//' @return NumericVector of n random variates.
//' @export
// [[Rcpp::export]]
NumericVector rsumgamma_series_cpp(int n, double sumalpha, double tau0, int i0,
                                  NumericVector cw_lower) {
  if (n <= 0) return NumericVector(0);
  
  const int kw = cw_lower.size() - 1;
  const double* p_cw = REAL(cw_lower);
  const double scale0 = 1.0 / tau0;
  NumericVector out(n);
  
  for (int idx = 0; idx < n; ++idx) {
    double u = R::unif_rand();
    // Binary search in cumulative weights [0, 1]
    const double* it = std::upper_bound(p_cw, p_cw + kw + 1, u);
    int j = static_cast<int>(it - p_cw) - 1;
    if (j < 0) j = 0;
    if (j >= kw) j = kw - 1;
    
    double ak = sumalpha + (i0 - 1) + j;
    out[idx] = R::rgamma(ak, scale0);
  }
  
  return out;
}

