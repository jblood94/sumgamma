#' Sum of Gammas Distribution Object
#'
#' Constructor for an S3 distribution object representing the sum of independent
#' Gamma random variables.
#'
#' @param shape Numeric vector of shape parameters.
#' @param rate Numeric vector of rate parameters.
#' @param scale Optional numeric vector of scale parameters (1 / rate).
#' @param method Character, underlying engine: \code{"auto"}, \code{"spline"},
#'   \code{"series"}, or \code{"integral"}.
#' @param ... Additional arguments passed to method constructor.
#' @return An S3 object of class \code{"sumgamma"}.
#' @export
sumgamma <- function(shape, rate, scale = 1 / rate,
                     method = c("auto", "spline", "series", "integral"), ...) {
  method <- match.arg(method)
  params <- prepare_gamma_params(shape, rate, scale)
  moms <- sumgamma_moments(shape = params$shape, rate = params$rate)
  
  actual_method <- if (method == "auto") "spline" else method
  
  closure <- switch(
    actual_method,
    spline = sumgamma_spline(shape = params$shape, rate = params$rate, ...),
    series = sumgamma_series(shape = params$shape, rate = params$rate, ...),
    integral = sumgamma_integral(shape = params$shape, rate = params$rate, ...)
  )
  
  obj <- list(
    params = params,
    moments = moms,
    method = actual_method,
    closure = closure,
    d = closure$d,
    p = closure$p,
    q = closure$q,
    r = closure$r
  )
  class(obj) <- "sumgamma"
  obj
}

#' @export
print.sumgamma <- function(x, ...) {
  cat(sprintf("Sum of %d Independent Gamma Random Variables\n", x$params$K))
  cat(sprintf("  Mean:     %0.4f\n", x$moments$mean))
  cat(sprintf("  Std Dev:  %0.4f\n", x$moments$sd))
  cat(sprintf("  Skewness: %0.4f\n", x$moments$skewness))
  cat(sprintf("  Engine:   %s\n", x$method))
  invisible(x)
}

#' @export
summary.sumgamma <- function(object, ...) {
  q_vals <- object$q(c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
  names(q_vals) <- c("1%", "5%", "25%", "50%", "75%", "95%", "99%")
  
  res <- list(
    components = object$params$K,
    sum_shape = sum(object$params$shape),
    moments = object$moments,
    quantiles = q_vals,
    method = object$method
  )
  class(res) <- "summary.sumgamma"
  res
}

#' @export
print.summary.sumgamma <- function(x, ...) {
  cat(sprintf("Summary of Sum of %d Gamma Distributions\n", x$components))
  cat("--------------------------------------------------\n")
  cat(sprintf("  Theoretical Mean:     %0.4f\n", x$moments$mean))
  cat(sprintf("  Theoretical Variance: %0.4f\n", x$moments$var))
  cat(sprintf("  Theoretical Std Dev:  %0.4f\n", x$moments$sd))
  cat(sprintf("  Skewness:             %0.4f\n", x$moments$skewness))
  cat(sprintf("  Excess Kurtosis:      %0.4f\n", x$moments$kurtosis_excess))
  cat("\nKey Quantiles:\n")
  print(round(x$quantiles, 4))
  invisible(x)
}

#' Plot Density and Cumulative Distribution for Sum of Gammas
#'
#' @param x An object of class \code{"sumgamma"}.
#' @param which Character or integer specifying \code{"density"}, \code{"cdf"}, or \code{"both"}.
#' @param n_points Integer, number of evaluation points along the curve.
#' @param ... Additional graphical parameters passed to plot.
#' @export
plot.sumgamma <- function(x, which = c("both", "density", "cdf"), n_points = 200L, ...) {
  which <- match.arg(which)
  
  mu <- x$moments$mean
  sd <- x$moments$sd
  
  # Bounds from ~0.1% to 99.9%
  q_lo <- max(1e-8, x$q(0.001))
  q_hi <- x$q(0.999)
  
  xs <- seq(q_lo, q_hi, length.out = n_points)
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "both") {
    graphics::par(mfrow = c(1, 2))
  }
  
  if (which %in% c("both", "density")) {
    ys_d <- x$d(xs)
    graphics::plot(xs, ys_d, type = "l", lwd = 2, col = "steelblue",
                   xlab = "x", ylab = "Density",
                   main = sprintf("Density: Sum of %d Gammas", x$params$K), ...)
    graphics::abline(v = mu, col = "firebrick", lty = 2, lwd = 1.5)
    graphics::grid()
  }
  
  if (which %in% c("both", "cdf")) {
    ys_p <- x$p(xs)
    graphics::plot(xs, ys_p, type = "l", lwd = 2, col = "darkorange",
                   xlab = "x", ylab = "Cumulative Probability",
                   main = sprintf("CDF: Sum of %d Gammas", x$params$K), ...)
    graphics::abline(v = mu, col = "firebrick", lty = 2, lwd = 1.5)
    graphics::abline(h = c(0, 1), col = "gray70", lty = 3)
    graphics::grid()
  }
  
  invisible(x)
}

