// =============================================================================
// File: calc_psumgamma.cpp
// Description:
//   Exact double-precision evaluation of the CDF / Survival Function and PDF
//   of a sum of K independent, non-identical Gamma random variables:
//     X = sum_{i=1}^K Gamma(shape = alpha_i, scale = b_i)
//
//   Mathematical Architecture:
//     1. Unified Saddlepoint Contour Integration:
//        Shifts the Bromwich / Gil-Pelaez inversion contour to s = t_hat + i*z,
//        where K'(t_hat) = sum_{i=1}^K alpha_i * b_i / (1 - b_i * t_hat) = x.
//        - For x > mu (t_hat > 0): Computes the upper tail S(x) directly.
//        - For x < mu (t_hat < 0): Crosses the s = 0 pole and computes the
//          lower tail F(x) = 1 - S(x) directly.
//        - For x == mu (t_hat == 0): Uses the exact Gil-Pelaez limit at t_hat = 0.
//     2. Singularity-Adapted Hyperbolic Quadrature on [0, z_0^{(0)}]:
//        Uses z = s_scale * sinh(w) with s_scale = min(|t_hat|, 1 / b_tilde_max)
//        to analytically resolve both the Lorentzian pole at z = +/- i*|t_hat|
//        and the nearest branch point at z = +/- i / b_tilde_max.
//     3. Geometric Midpoint Splitting + Wynn's Epsilon Algorithm on Tail Cycles:
//        Integrates between exact roots of the strictly monotonic saddlepoint
//        phase, splitting wide initial segments at sqrt(z_curr * z_next) and
//        accelerating the alternating series via Wynn's epsilon algorithm.
//     4. O(1) Suffix-Moment Taylor Acceleration:
//        Precomputes suffix moments of alpha_i * b_tilde_i^m (m = 2..11) once
//        per x so that all components with b_tilde_i * z <= 0.025 are evaluated
//        in O(1) time without calling std::log1p or std::atan.
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include "gauss_legendre.h"

using namespace Rcpp;

// -----------------------------------------------------------------------------
// SuffixMoments:
//   Precomputes suffix sums M[m][i] = sum_{j=i}^{K-1} alpha_j * scale_j^m
//   for m = 2 .. 11 when `scale` is sorted in descending order.
// -----------------------------------------------------------------------------
struct FastKernelEvaluator {
  size_t K;
  const double* alpha;
  const double* a2;
  const double* scale;
  bool is_descending;
  std::array<std::vector<double>, 12> M;
  
  FastKernelEvaluator(size_t K_, const std::vector<double>& alpha_vec,
                      const std::vector<double>& a2_vec,
                      const std::vector<double>& scale_vec)
    : K(K_), alpha(alpha_vec.data()), a2(a2_vec.data()), scale(scale_vec.data()) {
    is_descending = std::is_sorted(scale_vec.begin(), scale_vec.end(), std::greater<double>());
    if (is_descending && K >= 16) {
      for (int m = 2; m <= 11; ++m) {
        M[m].resize(K + 1, 0.0);
      }
      for (size_t i = K; i-- > 0;) {
        double s = scale[i];
        double s2 = s * s;
        double term = alpha[i] * s2; // m = 2
        M[2][i] = M[2][i + 1] + term;
        for (int m = 3; m <= 11; ++m) {
          term *= s;
          M[m][i] = M[m][i + 1] + term;
        }
      }
    }
  }
  
  // Finds first index i_cut where scale[i_cut] * z <= 0.025 via binary search
  inline size_t find_cutoff(double z) const {
    if (!is_descending || K < 16 || z <= 0.0) return K;
    const double thresh = 0.025 / z;
    if (scale[0] <= thresh) return 0;
    if (scale[K - 1] > thresh) return K;
    size_t lo = 0, hi = K - 1;
    while (lo < hi) {
      size_t mid = (lo + hi) >> 1;
      if (scale[mid] <= thresh) hi = mid;
      else lo = mid + 1;
    }
    return lo;
  }
  
  // Evaluates base phase: sum_{i=0}^{K-1} alpha_i * (scale_i * z - atan(scale_i * z))
  inline double eval_phi_base(double z) const {
    size_t i_cut = find_cutoff(z);
    double phi = 0.0;
    for (size_t i = 0; i < i_cut; ++i) {
      double btz = scale[i] * z;
      if (btz < 1e-3) {
        double btz2 = btz * btz;
        phi += alpha[i] * (btz2 * btz) * ((1.0 / 3.0) - btz2 * (0.2 - btz2 * (1.0 / 7.0)));
      } else {
        phi += alpha[i] * (btz - std::atan(btz));
      }
    }
    if (i_cut < K) {
      double z2 = z * z;
      double z4 = z2 * z2;
      double z6 = z4 * z2;
      double z8 = z4 * z4;
      double z10 = z8 * z2;
      phi += z * ((1.0 / 3.0) * M[3][i_cut] * z2
                    - (1.0 / 5.0) * M[5][i_cut] * z4
                    + (1.0 / 7.0) * M[7][i_cut] * z6
                    - (1.0 / 9.0) * M[9][i_cut] * z8
                    + (1.0 / 11.0) * M[11][i_cut] * z10);
    }
    return phi;
  }
  
  // Evaluates both phi_base(z) and its derivative d(phi_base)/dz simultaneously for Newton-Raphson
  inline void eval_phi_and_df_base(double z, double& phi, double& df) const {
    size_t i_cut = find_cutoff(z);
    phi = 0.0;
    df = 0.0;
    for (size_t i = 0; i < i_cut; ++i) {
      double btz = scale[i] * z;
      double btz2 = btz * btz;
      if (btz < 1e-3) {
        double btz2_sub = btz2 * (0.2 - btz2 * (1.0 / 7.0));
        phi += alpha[i] * (btz2 * btz) * ((1.0 / 3.0) - btz2_sub);
      } else {
        phi += alpha[i] * (btz - std::atan(btz));
      }
      df += alpha[i] * scale[i] * btz2 / (1.0 + btz2);
    }
    if (i_cut < K) {
      double z2 = z * z;
      double z4 = z2 * z2;
      double z6 = z4 * z2;
      double z8 = z4 * z4;
      double z10 = z8 * z2;
      phi += z * ((1.0 / 3.0) * M[3][i_cut] * z2
                    - (1.0 / 5.0) * M[5][i_cut] * z4
                    + (1.0 / 7.0) * M[7][i_cut] * z6
                    - (1.0 / 9.0) * M[9][i_cut] * z8
                    + (1.0 / 11.0) * M[11][i_cut] * z10);
      df += M[3][i_cut] * z2 - M[5][i_cut] * z4 + M[7][i_cut] * z6
            - M[9][i_cut] * z8 + M[11][i_cut] * z10;
    }
  }
  
  // Fused evaluation of u(z) = sum a2_i * log1p((scale_i*z)^2) and phi_base(z)
  inline void eval_u_and_phi_base(double z, double& u_z, double& phi) const {
    size_t i_cut = find_cutoff(z);
    u_z = 0.0;
    phi = 0.0;
    for (size_t i = 0; i < i_cut; ++i) {
      double btz = scale[i] * z;
      double btz2 = btz * btz;
      u_z += a2[i] * std::log1p(btz2);
      if (btz < 1e-3) {
        phi += alpha[i] * (btz2 * btz) * ((1.0 / 3.0) - btz2 * (0.2 - btz2 * (1.0 / 7.0)));
      } else {
        phi += alpha[i] * (btz - std::atan(btz));
      }
    }
    if (i_cut < K) {
      double z2 = z * z;
      double z4 = z2 * z2;
      double z6 = z4 * z2;
      double z8 = z4 * z4;
      double z10 = z8 * z2;
      u_z += 0.5 * (M[2][i_cut] * z2
                      - 0.5 * M[4][i_cut] * z4
                      + (1.0 / 3.0) * M[6][i_cut] * z6
                      - 0.25 * M[8][i_cut] * z8
                      + 0.2 * M[10][i_cut] * z10);
      phi += z * ((1.0 / 3.0) * M[3][i_cut] * z2
                    - (1.0 / 5.0) * M[5][i_cut] * z4
                    + (1.0 / 7.0) * M[7][i_cut] * z6
                    - (1.0 / 9.0) * M[9][i_cut] * z8
                    + (1.0 / 11.0) * M[11][i_cut] * z10);
    }
  }
};

// -----------------------------------------------------------------------------
// find_saddlepoint_delta_right:
//   Solves K'(t_hat) = x for x > mu in terms of the gap variable
//   delta = 1 - b_max * t_hat in (0, 1) to prevent 1 - b_max * t_hat underflow.
// -----------------------------------------------------------------------------
double find_saddlepoint_delta_right(double x, double b_max,
                                    const std::vector<double>& alpha,
                                    const std::vector<double>& r_vec,
                                    const std::vector<double>& ab) {
  const size_t K = alpha.size();
  const double* r = r_vec.data();
  const double* p_ab = ab.data();
  
  double alpha_max_sum = 0.0;
  for (size_t i = 0; i < K; ++i) {
    if (r[i] >= 1.0 - 1e-14) alpha_max_sum += alpha[i];
  }
  
  double lo = std::max(1e-300, 0.5 * (alpha_max_sum * b_max) / x);
  double hi = 1.0;
  
  while (true) {
    double kp = 0.0;
    for (size_t i = 0; i < K; ++i) {
      kp += p_ab[i] / ((1.0 - r[i]) + r[i] * lo);
    }
    if (kp >= x || lo <= 1e-300) break;
    lo *= 0.5;
  }
  
  double delta = std::sqrt(lo * hi);
  for (int iter = 0; iter < 80; ++iter) {
    double f = -x;
    double df = 0.0;
    for (size_t i = 0; i < K; ++i) {
      double inv = 1.0 / ((1.0 - r[i]) + r[i] * delta);
      double term = p_ab[i] * inv;
      f += term;
      df -= term * r[i] * inv;
    }
    if (std::abs(f) <= 1e-14 * x) break;
    
    if (f > 0.0) lo = delta;
    else hi = delta;
    
    double next_delta = delta - f / df;
    if (next_delta <= lo || next_delta >= hi) {
      next_delta = std::sqrt(lo * hi);
    }
    if (std::abs(next_delta - delta) <= 1e-15 * delta) {
      delta = next_delta;
      break;
    }
    delta = next_delta;
  }
  return delta;
}

// -----------------------------------------------------------------------------
// find_saddlepoint_abs_t_left:
//   Solves K'(-abs_t) = x for x < mu where abs_t = -t_hat > 0.
// -----------------------------------------------------------------------------
double find_saddlepoint_abs_t_left(double x, double b_max,
                                   const std::vector<double>& b,
                                   const std::vector<double>& ab) {
  const size_t K = b.size();
  const double* p_b = b.data();
  const double* p_ab = ab.data();
  
  double lo = 0.0;
  double hi = 1.0 / b_max;
  
  while (true) {
    double kp = 0.0;
    for (size_t i = 0; i < K; ++i) {
      kp += p_ab[i] / (1.0 + p_b[i] * hi);
    }
    if (kp <= x) break;
    lo = hi;
    hi *= 2.0;
  }
  
  double abs_t = 0.5 * (lo + hi);
  for (int iter = 0; iter < 80; ++iter) {
    double f = -x;
    double df = 0.0;
    for (size_t i = 0; i < K; ++i) {
      double inv = 1.0 / (1.0 + p_b[i] * abs_t);
      double term = p_ab[i] * inv;
      f += term;
      df -= term * p_b[i] * inv;
    }
    if (std::abs(f) <= 1e-14 * x) break;
    
    if (f > 0.0) lo = abs_t;
    else hi = abs_t;
    
    double next_t = abs_t - f / df;
    if (next_t <= lo || next_t >= hi) {
      next_t = (lo > 0.0) ? std::sqrt(lo * hi) : 0.5 * (lo + hi);
    }
    if (std::abs(next_t - abs_t) <= 1e-15 * abs_t) {
      abs_t = next_t;
      break;
    }
    abs_t = next_t;
  }
  return abs_t;
}

// -----------------------------------------------------------------------------
// accelerate_wynn_epsilon:
//   Computes the higher-order Shanks transformation of partial sums `II`
//   via Wynn's epsilon algorithm using a stack-allocated array.
// -----------------------------------------------------------------------------
double accelerate_wynn_epsilon(const std::vector<double>& II) {
  const size_t n = II.size();
  if (n == 0) return 0.0;
  if (n < 3) return II.back();
  
  constexpr size_t MAX_N = 64;
  double e[MAX_N][MAX_N] = {{0.0}};
  for (size_t i = 0; i < n; ++i) e[i][0] = II[i];
  
  double best_val = II.back();
  double min_score = std::abs(II[n - 1] - II[n - 2]);
  
  for (size_t k = 1; k < n; ++k) {
    for (size_t i = 0; i < n - k; ++i) {
      double diff = e[i + 1][k - 1] - e[i][k - 1];
      if (std::abs(diff) < 1e-30 || !std::isfinite(diff)) {
        e[i][k] = 1e300;
      } else {
        double prev = (k >= 2) ? e[i + 1][k - 2] : 0.0;
        e[i][k] = prev + 1.0 / diff;
      }
    }
    if (k % 2 == 0 && k <= 14) {
      for (size_t i = 1; i < n - k; ++i) {
        if (std::abs(e[i][k]) < 1e100 && std::abs(e[i - 1][k]) < 1e100) {
          double d_vert = std::abs(e[i][k] - e[i - 1][k]);
          double d_diag = std::abs(e[i][k] - e[i + 1][k - 2]);
          double score = d_vert + d_diag;
          if (score <= min_score) {
            min_score = score;
            best_val = e[i][k];
          }
        }
      }
    }
  }
  return best_val;
}

//' Evaluate Gamma Sum Distribution via Saddlepoint Contour Integration
//'
//' Evaluates the survival function or CDF of a sum of independent Gammas
//' via saddlepoint-shifted contour integration.
//'
//' @param x_in NumericVector of evaluation points.
//' @param mu Numeric, theoretical mean sum(alpha / tau).
//' @param nn Integer, number of Gauss-Legendre quadrature points per segment.
//' @param alpha_in NumericVector, shape parameters.
//' @param b_in NumericVector, scale parameters (1 / tau) sorted descending.
//' @param ab_in NumericVector, alpha * b.
//' @param a2_in NumericVector, alpha * 0.5.
//' @param nodes_x NumericVector, quadrature nodes in (0, 1). If empty, computed automatically.
//' @param nodes_w NumericVector, quadrature weights. If empty, computed automatically.
//' @param log_p Logical, if TRUE, returns log probabilities.
//' @param lower_tail Logical, if TRUE, returns P(X <= x); if FALSE, returns P(X > x).
//' @return NumericVector of probabilities.
//' @export
// [[Rcpp::export]]
NumericVector psumgamma1_cpp(NumericVector x_in, double mu, int nn,
                             NumericVector alpha_in, NumericVector b_in,
                             NumericVector ab_in, NumericVector a2_in,
                             NumericVector nodes_x, NumericVector nodes_w,
                             bool log_p = true, bool lower_tail = false) {
  
  std::vector<double> alpha = as<std::vector<double>>(alpha_in);
  std::vector<double> b = as<std::vector<double>>(b_in);
  std::vector<double> ab = as<std::vector<double>>(ab_in);
  std::vector<double> a2 = as<std::vector<double>>(a2_in);
  const size_t K = alpha.size();
  
  double b_max = 0.0;
  double sum_alpha = 0.0;
  for (double ai : alpha) sum_alpha += ai;
  const int max_ni = (sum_alpha < 1.0) ? 60 : 35;
  for (double bi : b) if (bi > b_max) b_max = bi;
  
  // Precompute r[i] = b[i] / b_max
  std::vector<double> r_vec(K);
  for (size_t i = 0; i < K; ++i) r_vec[i] = b[i] / b_max;
  
  // Auto-generate quadrature nodes if not provided
  std::vector<double> auto_nodes_x, auto_nodes_w;
  const double* p_nodes_x = nullptr;
  const double* p_nodes_w = nullptr;
  
  if (nodes_x.size() == nn && nodes_w.size() == nn) {
    p_nodes_x = REAL(nodes_x);
    p_nodes_w = REAL(nodes_w);
  } else {
    sumgamma::compute_gauss_legendre_01(nn, auto_nodes_x, auto_nodes_w);
    p_nodes_x = auto_nodes_x.data();
    p_nodes_w = auto_nodes_w.data();
  }
  
  const int n_out = x_in.size();
  NumericVector out(n_out);
  std::vector<double> b_tilde(K);
  std::vector<double> II;
  II.reserve(max_ni);
  
  for (int idx = 0; idx < n_out; ++idx) {
    double x = x_in[idx];
    if (x <= 0.0) {
      double p_val = lower_tail ? 0.0 : 1.0;
      out[idx] = log_p ? (p_val > 0.0 ? 0.0 : -R_PosInf) : p_val;
      continue;
    }
    if (std::isinf(x)) {
      double p_val = lower_tail ? 1.0 : 0.0;
      out[idx] = log_p ? (p_val > 0.0 ? 0.0 : -R_PosInf) : p_val;
      continue;
    }
    
    // Case 0: At the mean
    if (std::abs(x - mu) <= 1e-12 * mu) {
      FastKernelEvaluator ker_mu(K, alpha, a2, b);
      
      auto find_root_mu = [&](int k, double z_low) {
        double target_y = (k + 1) * M_PI;
        double lo = z_low;
        double hi = (k == 0) ? std::max(1.0 / mu, 1.0 / b_max) : z_low + (M_PI / mu);
        while (ker_mu.eval_phi_base(hi) < target_y) {
          lo = hi;
          hi *= 2.0;
        }
        double z = 0.5 * (lo + hi);
        for (int iter = 0; iter < 80; ++iter) {
          double phi, df;
          ker_mu.eval_phi_and_df_base(z, phi, df);
          double f = phi - target_y;
          if (std::abs(f) < 1e-14) break;
          if (f < 0.0) lo = z;
          else hi = z;
          double next_z = (df > 1e-15) ? (z - f / df) : 0.5 * (lo + hi);
          if (next_z <= lo || next_z >= hi) next_z = 0.5 * (lo + hi);
          if (std::abs(next_z - z) <= 1e-15 * z) { z = next_z; break; }
          z = next_z;
        }
        return z;
      };
      
      double z_curr = find_root_mu(0, 0.0);
      double s_scale = 1.0 / b_max;
      double w_max = std::asinh(z_curr / s_scale);
      double seg0_int = 0.0;
      for (int p = 0; p < 2; ++p) {
        double w_start = 0.5 * p * w_max;
        double dw = 0.5 * w_max;
        double panel_sum = 0.0;
        for (int j = 0; j < nn; ++j) {
          double w = w_start + dw * p_nodes_x[j];
          double z = s_scale * std::sinh(w);
          double dz_dw = s_scale * std::cosh(w);
          double u_z, phi;
          ker_mu.eval_u_and_phi_base(z, u_z, phi);
          panel_sum += (std::exp(-u_z) * (std::sin(phi) / z)) * dz_dw * p_nodes_w[j];
        }
        seg0_int += 0.5 * dw * panel_sum;
      }
      
      II.clear();
      double cumsum = 0.0;
      for (int k = 0; k < max_ni - 1; ++k) {
        double z_next = find_root_mu(k + 1, z_curr);
        double dz = z_next - z_curr;
        double seg_int = 0.0;
        for (int j = 0; j < nn; ++j) {
          double z = z_curr + dz * p_nodes_x[j];
          double u_z, phi;
          ker_mu.eval_u_and_phi_base(z, u_z, phi);
          seg_int += (std::exp(-u_z) * (std::sin(phi) / z)) * p_nodes_w[j];
        }
        seg_int *= 0.5 * dz;
        cumsum += seg_int;
        II.push_back(cumsum);
        z_curr = z_next;
        if (std::abs(seg_int) < 1e-16 && k >= 2) break;
      }
      double S_val = 0.5 - (seg0_int + accelerate_wynn_epsilon(II)) / M_PI;
      double P_val = lower_tail ? (1.0 - S_val) : S_val;
      out[idx] = log_p ? std::log(P_val) : P_val;
      continue;
    }
    
    // Case 1 & 2: Unified Saddlepoint Contour (x > mu or x < mu)
    const bool is_right = (x > mu);
    double abs_t = 0.0;
    double psi = 0.0;
    double b_tilde_max = 0.0;
    
    if (is_right) {
      double delta = find_saddlepoint_delta_right(x, b_max, alpha, r_vec, ab);
      abs_t = (1.0 - delta) / b_max;
      psi = -abs_t * x;
      for (size_t i = 0; i < K; ++i) {
        double one_minus_bt = (1.0 - r_vec[i]) + r_vec[i] * delta;
        psi -= alpha[i] * std::log(one_minus_bt);
        b_tilde[i] = b[i] / one_minus_bt;
        if (b_tilde[i] > b_tilde_max) b_tilde_max = b_tilde[i];
      }
    } else {
      abs_t = find_saddlepoint_abs_t_left(x, b_max, b, ab);
      psi = abs_t * x;
      for (size_t i = 0; i < K; ++i) {
        double bt = b[i] * abs_t;
        psi -= alpha[i] * std::log1p(bt);
        b_tilde[i] = b[i] / (1.0 + bt);
        if (b_tilde[i] > b_tilde_max) b_tilde_max = b_tilde[i];
      }
    }
    
    FastKernelEvaluator ker(K, alpha, a2, b_tilde);
    
    auto eval_total_phase = [&](double z) {
      double base = ker.eval_phi_base(z);
      double angle = std::atan2(z, abs_t);
      return is_right ? (base + angle) : (base - angle);
    };
    
    auto find_root_sp = [&](int k, double z_low) {
      double target_y = (k + 0.5) * M_PI;
      double lo = z_low;
      double hi = (k == 0) ? std::max(abs_t, 1.0 / x) : z_low + (M_PI / x);
      while (eval_total_phase(hi) < target_y) {
        lo = hi;
        hi *= 2.0;
      }
      double z = 0.5 * (lo + hi);
      for (int iter = 0; iter < 80; ++iter) {
        double base_phi, base_df;
        ker.eval_phi_and_df_base(z, base_phi, base_df);
        double angle = std::atan2(z, abs_t);
        double f = (is_right ? (base_phi + angle) : (base_phi - angle)) - target_y;
        if (std::abs(f) < 1e-14) break;
        if (f < 0.0) lo = z;
        else hi = z;
        
        double df = base_df + (is_right ? 1.0 : -1.0) * (abs_t / (abs_t * abs_t + z * z));
        double next_z = (df > 1e-15) ? (z - f / df) : 0.5 * (lo + hi);
        if (next_z <= lo || next_z >= hi) next_z = 0.5 * (lo + hi);
        if (std::abs(next_z - z) <= 1e-15 * z) {
          z = next_z;
          break;
        }
        z = next_z;
      }
      return z;
    };
    
    double z_curr = find_root_sp(0, 0.0);
    double s_scale = std::min(abs_t, 1.0 / b_tilde_max);
    double w_max = std::asinh(z_curr / s_scale);
    double seg0_int = 0.0;
    
    for (int p = 0; p < 2; ++p) {
      double w_start = 0.5 * p * w_max;
      double dw = 0.5 * w_max;
      double panel_sum = 0.0;
      for (int j = 0; j < nn; ++j) {
        double w = w_start + dw * p_nodes_x[j];
        double z = s_scale * std::sinh(w);
        double dz_dw = s_scale * std::cosh(w);
        
        double u_z, base_phi;
        ker.eval_u_and_phi_base(z, u_z, base_phi);
        double angle = std::atan2(z, abs_t);
        double phi_t = is_right ? (base_phi + angle) : (base_phi - angle);
        double kernel = std::cos(phi_t) / std::hypot(abs_t, z);
        panel_sum += std::exp(-u_z) * kernel * dz_dw * p_nodes_w[j];
      }
      seg0_int += 0.5 * dw * panel_sum;
    }
    
    auto integrate_linear_panel = [&](double za, double zb) {
      double dz = zb - za;
      double s = 0.0;
      for (int j = 0; j < nn; ++j) {
        double z = za + dz * p_nodes_x[j];
        double u_z, base_phi;
        ker.eval_u_and_phi_base(z, u_z, base_phi);
        double angle = std::atan2(z, abs_t);
        double phi_t = is_right ? (base_phi + angle) : (base_phi - angle);
        double kernel = std::cos(phi_t) / std::hypot(abs_t, z);
        s += std::exp(-u_z) * kernel * p_nodes_w[j];
      }
      return 0.5 * dz * s;
    };
    
    II.clear();
    double cumsum = 0.0;
    for (int k = 0; k < max_ni - 1; ++k) {
      double z_next = find_root_sp(k + 1, z_curr);
      double seg_int = 0.0;
      if (z_next > 2.0 * z_curr) {
        double z_mid = std::sqrt(z_curr * z_next);
        seg_int = integrate_linear_panel(z_curr, z_mid) + integrate_linear_panel(z_mid, z_next);
      } else {
        seg_int = integrate_linear_panel(z_curr, z_next);
      }
      cumsum += seg_int;
      II.push_back(cumsum);
      z_curr = z_next;
      if (std::abs(seg_int) < 1e-16 * seg0_int && k >= 2) break;
    }
    
    double total_J = seg0_int + accelerate_wynn_epsilon(II);
    double log_tail = psi + std::log(std::max(1e-300, total_J / M_PI));
    
    if (is_right != lower_tail) {
      out[idx] = log_p ? log_tail : std::exp(log_tail);
    } else {
      if (log_p) {
        out[idx] = (log_tail < -1e-16) ? std::log(-std::expm1(log_tail)) : -R_PosInf;
      } else {
        out[idx] = -std::expm1(log_tail);
      }
    }
  }
  
  return out;
}

//' Evaluate Gamma Sum Density via Saddlepoint Contour Integration
//'
//' Evaluates the probability density function (PDF) of a sum of independent Gammas
//' via saddlepoint-shifted contour integration of the inversion formula.
//'
//' @param x_in NumericVector of evaluation points.
//' @param mu Numeric, theoretical mean sum(alpha / tau).
//' @param nn Integer, number of Gauss-Legendre quadrature points per segment.
//' @param alpha_in NumericVector, shape parameters.
//' @param b_in NumericVector, scale parameters (1 / tau) sorted descending.
//' @param ab_in NumericVector, alpha * b.
//' @param a2_in NumericVector, alpha * 0.5.
//' @param nodes_x NumericVector, quadrature nodes in (0, 1). If empty, computed automatically.
//' @param nodes_w NumericVector, quadrature weights. If empty, computed automatically.
//' @param log_d Logical, if TRUE, returns log density.
//' @return NumericVector of densities.
//' @export
// [[Rcpp::export]]
NumericVector dsumgamma1_cpp(NumericVector x_in, double mu, int nn,
                             NumericVector alpha_in, NumericVector b_in,
                             NumericVector ab_in, NumericVector a2_in,
                             NumericVector nodes_x, NumericVector nodes_w,
                             bool log_d = false) {
  
  std::vector<double> alpha = as<std::vector<double>>(alpha_in);
  std::vector<double> b = as<std::vector<double>>(b_in);
  std::vector<double> ab = as<std::vector<double>>(ab_in);
  std::vector<double> a2 = as<std::vector<double>>(a2_in);
  const size_t K = alpha.size();
  
  double b_max = 0.0;
  double sum_alpha = 0.0;
  for (double ai : alpha) sum_alpha += ai;
  const int max_ni = (sum_alpha < 1.0) ? 60 : 35;
  for (double bi : b) if (bi > b_max) b_max = bi;
  
  std::vector<double> r_vec(K);
  for (size_t i = 0; i < K; ++i) r_vec[i] = b[i] / b_max;
  
  std::vector<double> auto_nodes_x, auto_nodes_w;
  const double* p_nodes_x = nullptr;
  const double* p_nodes_w = nullptr;
  
  if (nodes_x.size() == nn && nodes_w.size() == nn) {
    p_nodes_x = REAL(nodes_x);
    p_nodes_w = REAL(nodes_w);
  } else {
    sumgamma::compute_gauss_legendre_01(nn, auto_nodes_x, auto_nodes_w);
    p_nodes_x = auto_nodes_x.data();
    p_nodes_w = auto_nodes_w.data();
  }
  
  const int n_out = x_in.size();
  NumericVector out(n_out);
  std::vector<double> b_tilde(K);
  std::vector<double> II;
  II.reserve(max_ni);
  
  for (int idx = 0; idx < n_out; ++idx) {
    double x = x_in[idx];
    if (x <= 0.0) {
      out[idx] = log_d ? -R_PosInf : 0.0;
      continue;
    }
    if (std::isinf(x)) {
      out[idx] = log_d ? -R_PosInf : 0.0;
      continue;
    }
    
    // Saddlepoint solution t_hat
    double abs_t = 0.0;
    double psi = 0.0;
    double b_tilde_max = 0.0;
    bool at_mean = (std::abs(x - mu) <= 1e-12 * mu);
    
    if (at_mean) {
      psi = 0.0;
      for (size_t i = 0; i < K; ++i) {
        b_tilde[i] = b[i];
        if (b_tilde[i] > b_tilde_max) b_tilde_max = b_tilde[i];
      }
    } else if (x > mu) {
      double delta = find_saddlepoint_delta_right(x, b_max, alpha, r_vec, ab);
      abs_t = (1.0 - delta) / b_max;
      psi = -abs_t * x;
      for (size_t i = 0; i < K; ++i) {
        double one_minus_bt = (1.0 - r_vec[i]) + r_vec[i] * delta;
        psi -= alpha[i] * std::log(one_minus_bt);
        b_tilde[i] = b[i] / one_minus_bt;
        if (b_tilde[i] > b_tilde_max) b_tilde_max = b_tilde[i];
      }
    } else {
      abs_t = find_saddlepoint_abs_t_left(x, b_max, b, ab);
      psi = abs_t * x;
      for (size_t i = 0; i < K; ++i) {
        double bt = b[i] * abs_t;
        psi -= alpha[i] * std::log1p(bt);
        b_tilde[i] = b[i] / (1.0 + bt);
        if (b_tilde[i] > b_tilde_max) b_tilde_max = b_tilde[i];
      }
    }
    
    FastKernelEvaluator ker(K, alpha, a2, b_tilde);
    
    auto find_root_pdf = [&](int k, double z_low) {
      double target_y = (k + 0.5) * M_PI;
      double lo = z_low;
      double hi = (k == 0) ? std::max(1.0 / b_tilde_max, 1.0 / x) : z_low + (M_PI / x);
      while (ker.eval_phi_base(hi) < target_y) {
        lo = hi;
        hi *= 2.0;
      }
      double z = 0.5 * (lo + hi);
      for (int iter = 0; iter < 80; ++iter) {
        double base_phi, base_df;
        ker.eval_phi_and_df_base(z, base_phi, base_df);
        double f = base_phi - target_y;
        if (std::abs(f) < 1e-14) break;
        if (f < 0.0) lo = z;
        else hi = z;
        double next_z = (base_df > 1e-15) ? (z - f / base_df) : 0.5 * (lo + hi);
        if (next_z <= lo || next_z >= hi) next_z = 0.5 * (lo + hi);
        if (std::abs(next_z - z) <= 1e-15 * z) { z = next_z; break; }
        z = next_z;
      }
      return z;
    };
    
    // Panel 0: [0, z_0]
    double z_curr = find_root_pdf(0, 0.0);
    double s_scale = 1.0 / b_tilde_max;
    double w_max = std::asinh(z_curr / s_scale);
    double seg0_int = 0.0;
    
    for (int p = 0; p < 2; ++p) {
      double w_start = 0.5 * p * w_max;
      double dw = 0.5 * w_max;
      double panel_sum = 0.0;
      for (int j = 0; j < nn; ++j) {
        double w = w_start + dw * p_nodes_x[j];
        double z = s_scale * std::sinh(w);
        double dz_dw = s_scale * std::cosh(w);
        double u_z, base_phi;
        ker.eval_u_and_phi_base(z, u_z, base_phi);
        panel_sum += std::exp(-u_z) * std::cos(base_phi) * dz_dw * p_nodes_w[j];
      }
      seg0_int += 0.5 * dw * panel_sum;
    }
    
    auto integrate_linear_panel = [&](double za, double zb) {
      double dz = zb - za;
      double s = 0.0;
      for (int j = 0; j < nn; ++j) {
        double z = za + dz * p_nodes_x[j];
        double u_z, base_phi;
        ker.eval_u_and_phi_base(z, u_z, base_phi);
        s += std::exp(-u_z) * std::cos(base_phi) * p_nodes_w[j];
      }
      return 0.5 * dz * s;
    };
    
    II.clear();
    double cumsum = 0.0;
    for (int k = 0; k < max_ni - 1; ++k) {
      double z_next = find_root_pdf(k + 1, z_curr);
      double seg_int = 0.0;
      if (z_next > 2.0 * z_curr) {
        double z_mid = std::sqrt(z_curr * z_next);
        seg_int = integrate_linear_panel(z_curr, z_mid) + integrate_linear_panel(z_mid, z_next);
      } else {
        seg_int = integrate_linear_panel(z_curr, z_next);
      }
      cumsum += seg_int;
      II.push_back(cumsum);
      z_curr = z_next;
      if (std::abs(seg_int) < 1e-16 * seg0_int && k >= 2) break;
    }
    
    double total_J = seg0_int + accelerate_wynn_epsilon(II);
    double log_dens = psi + std::log(std::max(1e-300, total_J / M_PI));
    out[idx] = log_d ? log_dens : std::exp(log_dens);
  }
  
  return out;
}

