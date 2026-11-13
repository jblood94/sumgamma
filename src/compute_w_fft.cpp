#include <Rcpp.h>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// fft_radix2_ld:
//   In-place iterative Cooley-Tukey Radix-2 FFT over extended precision `long double`.
//   Uses precomputed half-circle twiddle table so forward and inverse transforms
//   share exact trigonometric evaluations without recurrence drift.
// -----------------------------------------------------------------------------
void fft_radix2_ld(std::vector<std::complex<long double>>& a,
                   const std::vector<std::complex<long double>>& twiddles,
                   bool invert) {
  const size_t n = a.size();
  
  // 1. Bit-reversal permutation
  for (size_t i = 1, j = 0; i < n; ++i) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  
  // 2. Butterfly stages using precomputed twiddle table
  for (size_t len = 2; len <= n; len <<= 1) {
    const size_t half = len >> 1;
    const size_t step = n / len;
    
    for (size_t i = 0; i < n; i += len) {
      for (size_t j = 0, k = 0; j < half; ++j, k += step) {
        std::complex<long double> w = twiddles[k];
        if (invert) w = std::conj(w);
        
        std::complex<long double> u = a[i + j];
        std::complex<long double> v = a[i + j + half] * w;
        a[i + j] = u + v;
        a[i + j + half] = u - v;
      }
    }
  }
  
  // 3. Rescale by 1 / n for inverse transform
  if (invert) {
    const long double inv_n = 1.0L / static_cast<long double>(n);
    for (std::complex<long double>& x : a) {
      x *= inv_n;
    }
  }
}

//' Compute Moschopoulos Infinite Mixture Weights via Radix-2 FFT
//'
//' Evaluates the first M mixture weights w_k of the Moschopoulos (1985) expansion
//' for the sum of independent Gamma random variables using an 80-bit
//' extended-precision Fast Fourier Transform.
//'
//' @param M Integer, power-of-2 FFT truncation length.
//' @param logv NumericVector, log(1 - tau / tau0) where tau is sorted ascending.
//' @param alpha NumericVector, shape parameters.
//' @param logtaun NumericVector, log(tau / tau0).
//' @param rel_tol Numeric, relative truncation tolerance for tail zeroing.
//' @return NumericVector of normalized mixture weights of length M.
//' @export
// [[Rcpp::export]]
NumericVector compute_w_fft_ld(int M, NumericVector logv, NumericVector alpha,
                               NumericVector logtaun, double rel_tol = 1e-17) {
  const int N = alpha.size();
  
  // Round FFT length up to the next power of 2
  size_t n_fft = 1;
  while (n_fft < static_cast<size_t>(M)) n_fft <<= 1;
  
  const double* p_logv = REAL(logv);
  const double* p_alpha = REAL(alpha);
  const double* p_logtaun = REAL(logtaun);
  
  std::vector<std::complex<long double>> buf(n_fft, std::complex<long double>(0.0L, 0.0L));
  long double g0 = 0.0L;
  
  std::vector<long double> v(N);
  std::vector<long double> current_val(N);
  
  for (int i = 0; i < N; ++i) {
    long double a_i = static_cast<long double>(p_alpha[i]);
    g0 += a_i * static_cast<long double>(p_logtaun[i]);
    v[i] = std::exp(static_cast<long double>(p_logv[i]));
    current_val[i] = a_i * v[i];
  }
  buf[0] = std::complex<long double>(g0, 0.0L);
  
  // Active-set shrinking: since tau is sorted ascending,
  // v[i] = 1 - tau[i]/tau0 is sorted descending (v[N-1] = 0).
  int active_N = N;
  while (active_N > 0 && current_val[active_N - 1] <= 0.0L) {
    --active_N;
  }
  
  for (size_t k = 1; k < n_fft && active_N > 0; ++k) {
    long double s_k = 0.0L;
    for (int i = 0; i < active_N; ++i) {
      s_k += current_val[i];
      current_val[i] *= v[i];
      if (current_val[i] < 1e-300L) current_val[i] = 0.0L;
    }
    buf[k] = std::complex<long double>(s_k / static_cast<long double>(k), 0.0L);
    
    // Shrink active_N from the right as faster-decaying components hit 0
    while (active_N > 0 && current_val[active_N - 1] <= 0.0L) {
      --active_N;
    }
  }
  
  // Precompute twiddle factors exp(-2*pi*i*j / n_fft)
  const size_t half_fft = n_fft >> 1;
  std::vector<std::complex<long double>> twiddles(half_fft);
  const long double pi_ld = std::acos(-1.0L);
  const long double ang = -2.0L * pi_ld / static_cast<long double>(n_fft);
  const long double s_half = std::sin(0.5L * ang);
  const std::complex<long double> w_step(-2.0L * s_half * s_half, std::sin(ang));
  
  std::complex<long double> w_cur(1.0L, 0.0L);
  for (size_t j = 0; j < half_fft; ++j) {
    if ((j & 255) == 0 && j > 0) {
      long double theta = ang * static_cast<long double>(j);
      w_cur = std::complex<long double>(std::cos(theta), std::sin(theta));
    }
    twiddles[j] = w_cur;
    w_cur += w_cur * w_step;
  }
  
  // Forward FFT: G = FFT(g)
  fft_radix2_ld(buf, twiddles, false);
  
  // Subtract max_real before exponentiating to prevent overflow
  long double max_real = buf[0].real();
  for (size_t k = 1; k < n_fft; ++k) {
    if (buf[k].real() > max_real) max_real = buf[k].real();
  }
  
  for (size_t k = 0; k < n_fft; ++k) {
    long double mag = std::exp(buf[k].real() - max_real);
    long double phase = buf[k].imag();
    buf[k] = std::complex<long double>(mag * std::cos(phase), mag * std::sin(phase));
  }
  
  // Inverse FFT: w = Re(IFFT(exp(G)))
  fft_radix2_ld(buf, twiddles, true);
  
  // Locate mode k_peak and peak weight w_max
  size_t k_peak = 0;
  long double w_max = 0.0L;
  for (int k = 0; k < M; ++k) {
    if (buf[k].real() > w_max) {
      w_max = buf[k].real();
      k_peak = static_cast<size_t>(k);
    }
  }
  
  // Zero out weights below rel_tol * w_max and clip right tail once it hits floor
  const long double thresh = static_cast<long double>(rel_tol) * w_max;
  long double sum_w = 0.0L;
  bool tail_cut = false;
  
  for (int k = 0; k < M; ++k) {
    long double wk = buf[k].real();
    if (static_cast<size_t>(k) > k_peak && wk < thresh) {
      tail_cut = true;
    }
    if (tail_cut || wk < thresh) {
      buf[k] = std::complex<long double>(0.0L, 0.0L);
    } else {
      sum_w += wk;
    }
  }
  
  NumericVector out(M);
  double* p_out = REAL(out);
  const long double inv_sum_w = (sum_w > 0.0L) ? (1.0L / sum_w) : 1.0L;
  for (int k = 0; k < M; ++k) {
    p_out[k] = static_cast<double>(buf[k].real() * inv_sum_w);
  }
  
  return out;
}

//' High-Level Moschopoulos Mixture Weights Generator
//'
//' Prepares and computes trimmed, normalized Moschopoulos mixture weights,
//' along with prefix and suffix cumulative sums for rapid tail evaluation.
//'
//' @param alpha NumericVector, shape parameters (grouped and matching tau).
//' @param tau NumericVector, rate parameters sorted in ascending order.
//' @param M Integer, power-of-2 FFT truncation size.
//' @param rel_tol Numeric, relative truncation tolerance.
//' @return A List containing \code{w}, \code{i0}, \code{sumalpha}, \code{tau0},
//' \code{cw_lower}, and \code{cw_upper}.
//' @keywords internal
//' @export
// [[Rcpp::export]]
List compute_mixture_weights_cpp(NumericVector alpha, NumericVector tau,
                                 int M, double rel_tol = 1e-17) {
  int N = alpha.size();
  if (N != tau.size() || N == 0) stop("alpha and tau must be non-empty vectors of equal length");
  
  double tau0 = tau[N - 1];
  NumericVector taun = tau / tau0;
  NumericVector logv(N);
  NumericVector logtaun(N);
  double sumalpha = 0.0;
  
  for (int i = 0; i < N; ++i) {
    sumalpha += alpha[i];
    logtaun[i] = std::log(taun[i]);
    logv[i] = (i == N - 1) ? -700.0 : std::log1p(-taun[i]);
  }
  
  NumericVector w_raw = compute_w_fft_ld(M, logv, alpha, logtaun, rel_tol);
  
  // Find first positive weight index i0 (1-based index)
  int i0 = 1;
  while (i0 <= M && w_raw[i0 - 1] <= 0.0) ++i0;
  
  // Find last positive weight index i1 (1-based index)
  int i1 = M;
  while (i1 >= i0 && w_raw[i1 - 1] <= 0.0) --i1;
  
  if (i0 > i1) {
    // Fallback if all weights zeroed (degenerate)
    i0 = 1; i1 = 1;
  }
  
  int kw = i1 - i0 + 1;
  NumericVector w(kw);
  double sum_trimmed = 0.0;
  for (int j = 0; j < kw; ++j) {
    w[j] = w_raw[i0 - 1 + j];
    sum_trimmed += w[j];
  }
  
  if (sum_trimmed > 0.0) {
    double inv_s = 1.0 / sum_trimmed;
    for (int j = 0; j < kw; ++j) w[j] *= inv_s;
  }
  
  // Cumulative sums:
  // cw_lower has length kw + 1: cw_lower[0] = 0, cw_lower[j] = sum_{m < j} w[m]
  // cw_upper has length kw + 1: cw_upper[j] = sum_{m >= j} w[m], cw_upper[kw] = 0
  NumericVector cw_lower(kw + 1);
  NumericVector cw_upper(kw + 1);
  
  double acc_low = 0.0;
  cw_lower[0] = 0.0;
  for (int j = 0; j < kw; ++j) {
    acc_low += w[j];
    cw_lower[j + 1] = acc_low;
  }
  
  double acc_up = 0.0;
  cw_upper[kw] = 0.0;
  for (int j = kw - 1; j >= 0; --j) {
    acc_up += w[j];
    cw_upper[j] = acc_up;
  }
  
  return List::create(
    Named("w") = w,
    Named("i0") = i0,
    Named("sumalpha") = sumalpha,
    Named("tau0") = tau0,
    Named("cw_lower") = cw_lower,
    Named("cw_upper") = cw_upper
  );
}

