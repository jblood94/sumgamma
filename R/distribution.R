#' Direct Monte Carlo Simulation for Sum of Gammas
#'
#' Generates n exact random variates by summing direct draws from each Gamma component.
#'
#' @param n Integer, number of samples to generate.
#' @param shape Numeric vector of shape parameters.
#' @param rate Numeric vector of rate parameters.
#' @param scale Optional numeric vector of scale parameters (1 / rate).
#' @return Numeric vector of n random variates.
#' @export
rsumgamma_direct <- function(n, shape, rate, scale = 1 / rate) {
  params <- prepare_gamma_params(shape, rate, scale)
  alpha <- params$shape
  tau <- params$rate
  K <- params$K
  
  if (n <= 0) return(numeric(0))
  
  # Check if tiny shapes present that risk underflowing standard rgamma
  if (any(alpha < 1e-4)) {
    # Log-space gamma generation trick
    x <- matrix(log(stats::rgamma(K * n, rep.int(alpha + 1, n))) -
                  stats::rexp(K * n, rep.int(alpha, n)) -
                  log(rep.int(tau, n)), K, n)
    m <- apply(x, 2, max)
    colSums(exp(sweep(x, 2, m, "-"))) * exp(m)
  } else {
    colSums(matrix(stats::rgamma(K * n, rep.int(alpha, n), rep.int(tau, n)), K, n))
  }
}

#' Probability Density Function for Sum of Gammas
#'
#' Evaluates the density function of the convolution of independent Gamma variables.
#'
#' @param x Numeric vector of evaluation points.
#' @param shape Numeric vector of shape parameters.
#' @param rate Numeric vector of rate parameters.
#' @param scale Optional numeric vector of scale parameters (1 / rate).
#' @param log Logical, if TRUE, returns log-densities.
#' @param method Character, evaluation method: \code{"auto"}, \code{"spline"},
#'   \code{"series"}, or \code{"integral"}.
#' @param ... Additional arguments passed to method constructor.
#' @return Numeric vector of densities.
#' @export
dsumgamma <- function(x, shape, rate, scale = 1 / rate, log = FALSE,
                      method = c("auto", "spline", "series", "integral"), ...) {
  method <- match.arg(method)
  if (method == "auto") {
    method <- if (length(x) > 50L) "spline" else "series"
  }
  
  factory <- switch(
    method,
    spline = sumgamma_spline(shape, rate, scale, ...),
    series = sumgamma_series(shape, rate, scale, ...),
    integral = sumgamma_integral(shape, rate, scale, ...)
  )
  factory$d(x, log = log)
}

#' Cumulative Distribution Function for Sum of Gammas
#'
#' Evaluates the cumulative distribution (or survival) function of a sum of Gammas.
#'
#' @param q Numeric vector of quantiles.
#' @param shape Numeric vector of shape parameters.
#' @param rate Numeric vector of rate parameters.
#' @param scale Optional numeric vector of scale parameters (1 / rate).
#' @param lower.tail Logical, if TRUE, returns P(X <= q); if FALSE, returns P(X > q).
#' @param log.p Logical, if TRUE, returns log probabilities.
#' @param method Character, evaluation method: \code{"auto"}, \code{"spline"},
#'   \code{"series"}, or \code{"integral"}.
#' @param ... Additional arguments passed to method constructor.
#' @return Numeric vector of probabilities.
#' @export
psumgamma <- function(q, shape, rate, scale = 1 / rate,
                      lower.tail = TRUE, log.p = FALSE,
                      method = c("auto", "spline", "series", "integral"), ...) {
  method <- match.arg(method)
  if (method == "auto") {
    method <- if (length(q) > 50L) "spline" else "series"
  }
  
  factory <- switch(
    method,
    spline = sumgamma_spline(shape, rate, scale, ...),
    series = sumgamma_series(shape, rate, scale, ...),
    integral = sumgamma_integral(shape, rate, scale, ...)
  )
  factory$p(q, lower.tail = lower.tail, log.p = log.p)
}

#' Quantile Function for Sum of Gammas
#'
#' Computes quantiles for the convolution of independent Gamma variables.
#'
#' @param p Numeric vector of probabilities.
#' @param shape Numeric vector of shape parameters.
#' @param rate Numeric vector of rate parameters.
#' @param scale Optional numeric vector of scale parameters (1 / rate).
#' @param lower.tail Logical, if TRUE, returns quantiles for P(X <= q); if FALSE, for P(X > q).
#' @param log.p Logical, if TRUE, probabilities are given as log(p).
#' @param method Character, evaluation method: \code{"auto"}, \code{"spline"}, or \code{"integral"}.
#' @param ... Additional arguments passed to method constructor.
#' @return Numeric vector of quantiles.
#' @export
qsumgamma <- function(p, shape, rate, scale = 1 / rate,
                      lower.tail = TRUE, log.p = FALSE,
                      method = c("auto", "spline", "integral"), ...) {
  method <- match.arg(method)
  if (method == "auto") method <- "spline"
  
  factory <- switch(
    method,
    spline = sumgamma_spline(shape, rate, scale, ...),
    integral = sumgamma_integral(shape, rate, scale, ...)
  )
  factory$q(p, lower.tail = lower.tail, log.p = log.p)
}

#' Random Variate Generation for Sum of Gammas
#'
#' Generates random variates from the convolution of independent Gamma variables.
#'
#' @param n Integer, number of observations.
#' @param shape Numeric vector of shape parameters.
#' @param rate Numeric vector of rate parameters.
#' @param scale Optional numeric vector of scale parameters (1 / rate).
#' @param method Character, sampling method: \code{"auto"}, \code{"spline"},
#'   \code{"series"}, or \code{"direct"}.
#' @param ... Additional arguments passed to method constructor.
#' @return Numeric vector of random variates.
#' @export
rsumgamma <- function(n, shape, rate, scale = 1 / rate,
                      method = c("auto", "spline", "series", "direct"), ...) {
  method <- match.arg(method)
  if (method == "auto") {
    params <- prepare_gamma_params(shape, rate, scale)
    method <- if (n <= 10000L && params$K <= 50L) "direct" else "spline"
  }
  
  if (method == "direct") {
    return(rsumgamma_direct(n, shape, rate, scale))
  }
  
  factory <- switch(
    method,
    spline = sumgamma_spline(shape, rate, scale, ...),
    series = sumgamma_series(shape, rate, scale, ...)
  )
  factory$r(n)
}

