#' Construct Moschopoulos Infinite Mixture Distribution Functions
#'
#' Factory function returning distribution closures (\code{d}, \code{p}, \code{q}, \code{r})
#' for the sum of independent, non-identical Gamma random variables evaluated via
#' 80-bit extended-precision FFT Moschopoulos series expansion.
#'
#' @param shape Numeric vector of positive shape parameters.
#' @param rate Numeric vector of positive rate parameters.
#' @param scale Optional numeric vector of positive scale parameters (\code{1 / rate}).
#' @param kmax Integer, maximum power-of-2 FFT truncation size (default \code{2^21} = 2,097,152).
#' @param rel_tol Numeric, relative truncation tolerance for mixture tail truncation (default \code{1e-17}).
#' @return A list of closures of class \code{c("sumgamma_series", "sumgamma_closure")}:
#' \describe{
#'   \item{\code{d(x, log = FALSE)}}{Probability density function.}
#'   \item{\code{p(q, lower.tail = TRUE, log.p = FALSE)}}{Cumulative distribution / survival function.}
#'   \item{\code{q(p, lower.tail = TRUE, log.p = FALSE)}}{Quantile function (via bracketed numerical inversion).}
#'   \item{\code{r(n)}}{Exact two-stage categorical mixture random variate generator.}
#' }
#' @export
sumgamma_series <- function(shape, rate, scale = 1 / rate,
                            kmax = 2097152L, rel_tol = 1e-17) {
  params <- prepare_gamma_params(shape, rate, scale)
  
  alpha <- params$shape
  tau <- params$rate
  mu <- params$mu
  sigma <- params$sigma
  sumalpha <- sum(alpha)
  
  # Ensure kmax is a power of 2
  M <- 2^ceiling(log2(max(1024L, as.integer(kmax))))
  
  mix <- compute_mixture_weights_cpp(alpha, tau, M, rel_tol = rel_tol)
  w <- mix$w
  cw_lower <- mix$cw_lower
  cw_upper <- mix$cw_upper
  i0 <- mix$i0
  tau0 <- mix$tau0
  
  d_fn <- function(x, log = FALSE) {
    dsumgamma_series_cpp(as.numeric(x), sumalpha, tau0, i0, w, log_d = as.logical(log))
  }
  
  p_fn <- function(q, lower.tail = TRUE, log.p = FALSE) {
    psumgamma_series_cpp(as.numeric(q), sumalpha, tau0, i0, w, cw_lower, cw_upper,
                         lower_tail = as.logical(lower.tail),
                         log_p = as.logical(log.p))
  }
  
  r_fn <- function(n) {
    rsumgamma_series_cpp(as.integer(n), sumalpha, tau0, i0, cw_lower)
  }
  
  q_fn <- function(p, lower.tail = TRUE, log.p = FALSE) {
    p_num <- as.numeric(p)
    res <- numeric(length(p_num))
    
    p_lin <- if (log.p) {
      if (lower.tail) exp(p_num) else -expm1(p_num)
    } else {
      if (lower.tail) p_num else 1.0 - p_num
    }
    
    for (i in seq_along(p_lin)) {
      pi_val <- p_lin[i]
      if (is.na(pi_val)) { res[i] <- NA_real_; next }
      if (pi_val <= 0.0) { res[i] <- 0.0; next }
      if (pi_val >= 1.0) { res[i] <- Inf; next }
      
      target_log <- if (log.p && lower.tail) p_num[i] else log(pi_val)
      
      lo <- max(1e-12 * mu, mu - 10 * sigma)
      while (p_fn(lo, lower.tail = TRUE, log.p = TRUE) > target_log && lo > 1e-200 * mu) {
        lo <- lo * 0.1
      }
      
      hi <- mu + 10 * sigma
      while (p_fn(hi, lower.tail = TRUE, log.p = TRUE) < target_log && hi < 1e200 * mu) {
        hi <- hi * 2.0
      }
      
      obj <- function(x) p_fn(x, lower.tail = TRUE, log.p = TRUE) - target_log
      res[i] <- stats::uniroot(obj, interval = c(lo, hi), tol = 1e-12)$root
    }
    res
  }
  
  out <- list(d = d_fn, p = p_fn, q = q_fn, r = r_fn)
  attr(out, "params") <- params
  attr(out, "engine") <- "series"
  attr(out, "mixture") <- list(
    w = w,
    i0 = i0,
    k = length(w),
    M = M,
    tau0 = tau0,
    sumalpha = sumalpha
  )
  class(out) <- c("sumgamma_series", "sumgamma_closure")
  out
}

