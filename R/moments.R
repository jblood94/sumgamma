#' Theoretical Moments and Cumulants for Sum of Gammas
#'
#' Computes the exact theoretical mean, variance, skewness, and excess kurtosis
#' in O(K) time using the additivity of cumulants for independent variables.
#'
#' @param shape Numeric vector of shape parameters.
#' @param rate Numeric vector of rate parameters.
#' @param scale Optional numeric vector of scale parameters (1 / rate).
#' @return A list containing:
#' \describe{
#'   \item{mean}{Theoretical mean: \eqn{\mu = \sum \alpha_i / \tau_i}.}
#'   \item{var}{Theoretical variance: \eqn{\sigma^2 = \sum \alpha_i / \tau_i^2}.}
#'   \item{sd}{Theoretical standard deviation: \eqn{\sigma}.}
#'   \item{skewness}{Standardized skewness: \eqn{\kappa_3 / \sigma^3}.}
#'   \item{kurtosis_excess}{Excess kurtosis: \eqn{\kappa_4 / \sigma^4}.}
#'   \item{cumulants}{First four cumulants \eqn{(\kappa_1, \kappa_2, \kappa_3, \kappa_4)}.}
#' }
#' @export
sumgamma_moments <- function(shape, rate, scale = 1 / rate) {
  params <- prepare_gamma_params(shape, rate, scale)
  b <- params$scale
  a <- params$shape
  
  k1 <- sum(a * b)
  k2 <- sum(a * b^2)
  k3 <- 2 * sum(a * b^3)
  k4 <- 6 * sum(a * b^4)
  
  sigma <- sqrt(k2)
  skew <- k3 / (sigma^3)
  kurt <- k4 / (sigma^4)
  
  list(
    mean = k1,
    var = k2,
    sd = sigma,
    skewness = skew,
    kurtosis_excess = kurt,
    cumulants = c(k1 = k1, k2 = k2, k3 = k3, k4 = k4)
  )
}

#' Characteristic Function for Sum of Gammas
#'
#' Evaluates the exact complex characteristic function of the convolution of
#' independent Gamma random variables:
#' \deqn{\phi_X(t) = \prod_{k=1}^K (1 - i t / \tau_k)^{-\alpha_k}}
#'
#' @param t Numeric vector of evaluation arguments.
#' @param shape Numeric vector of shape parameters.
#' @param rate Numeric vector of rate parameters.
#' @param scale Optional numeric vector of scale parameters (1 / rate).
#' @return Complex vector of characteristic function values.
#' @export
cf_sumgamma <- function(t, shape, rate, scale = 1 / rate) {
  params <- prepare_gamma_params(shape, rate, scale)
  a <- params$shape
  b <- params$scale
  
  res <- complex(length(t))
  for (j in seq_along(t)) {
    tj <- t[j]
    terms <- (complex(real = 1, imaginary = -b * tj))^(-a)
    res[j] <- prod(terms)
  }
  res
}

#' Diagnostic and Algorithm Recommender for Sum of Gammas
#'
#' Analyzes the parameter configuration to determine rate dispersion, tail decay,
#' and optimal computational architecture.
#'
#' @param shape Numeric vector of shape parameters.
#' @param rate Numeric vector of rate parameters.
#' @param scale Optional numeric vector of scale parameters (1 / rate).
#' @param kmax Integer, crossover threshold.
#' @return A list containing diagnostic metrics and recommended algorithms.
#' @export
sumgamma_diagnose <- function(shape, rate, scale = 1 / rate, kmax = 2097152L) {
  params <- prepare_gamma_params(shape, rate, scale)
  moms <- sumgamma_moments(shape = params$shape, rate = params$rate)
  
  K <- params$K
  rate_ratio <- max(params$rate) / min(params$rate)
  tau0 <- max(params$rate)
  sumalpha <- sum(params$shape)
  
  # Approximate Moschopoulos truncation size k
  test_val <- (sumalpha + kmax) / tau0
  gq <- gauss_legendre_nodes(20L)
  b <- params$scale
  a <- params$shape
  ab <- a * b
  a2 <- a * 0.5
  
  tail_log_p <- psumgamma1_cpp(test_val, params$mu, 20L,
                               a, b, ab, a2, gq$nodes, gq$weights,
                               log_p = TRUE, lower_tail = FALSE)
  
  can_use_series <- (tail_log_p < -35)
  rec_engine <- if (can_use_series) "series (FFT Mixture)" else "integral (Saddlepoint Inversion)"
  
  list(
    components = K,
    sum_shape = sumalpha,
    rate_range = c(min = min(params$rate), max = max(params$rate)),
    rate_disparity_ratio = rate_ratio,
    mean = moms$mean,
    sd = moms$sd,
    skewness = moms$skewness,
    kmax_tail_log_prob = tail_log_p,
    recommended_engine = rec_engine,
    recommended_general_factory = "sumgamma_spline"
  )
}

