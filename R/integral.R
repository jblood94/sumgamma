#' Construct Saddlepoint Contour Integral Distribution Functions
#'
#' Factory function returning distribution closures (\code{d}, \code{p}, \code{q})
#' for the sum of independent, non-identical Gamma random variables evaluated via
#' high-precision saddlepoint contour integration.
#'
#' @param shape Numeric vector of positive shape parameters.
#' @param rate Numeric vector of positive rate parameters.
#' @param scale Optional numeric vector of positive scale parameters (\code{1 / rate}).
#' @param nn Integer, number of Gauss-Legendre quadrature points per contour segment (default 40).
#' @return A list of closures of class \code{c("sumgamma_integral", "sumgamma_closure")}:
#' \describe{
#'   \item{\code{d(x, log = FALSE)}}{Probability density function.}
#'   \item{\code{p(q, lower.tail = TRUE, log.p = FALSE)}}{Cumulative distribution / survival function.}
#'   \item{\code{q(p, lower.tail = TRUE, log.p = FALSE)}}{Quantile function (via bracketed numerical inversion).}
#' }
#' @export
sumgamma_integral <- function(shape, rate, scale = 1 / rate, nn = 40L) {
  params <- prepare_gamma_params(shape, rate, scale)
  
  # Parameters with scale sorted descending
  # Since rate is sorted ascending, scale = 1/rate is already sorted descending
  b <- params$scale
  alpha <- params$shape
  ab <- alpha * b
  a2 <- alpha * 0.5
  mu <- params$mu
  sigma <- params$sigma
  
  gq <- gauss_legendre_nodes(as.integer(nn))
  nodes_x <- gq$nodes
  nodes_w <- gq$weights
  
  d_fn <- function(x, log = FALSE) {
    dsumgamma1_cpp(as.numeric(x), mu, as.integer(nn),
                   alpha, b, ab, a2, nodes_x, nodes_w,
                   log_d = as.logical(log))
  }
  
  p_fn <- function(q, lower.tail = TRUE, log.p = FALSE) {
    psumgamma1_cpp(as.numeric(q), mu, as.integer(nn),
                   alpha, b, ab, a2, nodes_x, nodes_w,
                   log_p = as.logical(log.p),
                   lower_tail = as.logical(lower.tail))
  }
  
  q_fn <- function(p, lower.tail = TRUE, log.p = FALSE) {
    p_num <- as.numeric(p)
    res <- numeric(length(p_num))
    
    # Standardize to lower linear probability
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
      
      # Bracket the root using conservative sigma bounds
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
  
  out <- list(d = d_fn, p = p_fn, q = q_fn)
  attr(out, "params") <- params
  attr(out, "engine") <- "integral"
  class(out) <- c("sumgamma_integral", "sumgamma_closure")
  out
}

