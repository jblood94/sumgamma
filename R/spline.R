#' Extract Strictly Monotonic Points for Hyman Spline Fitting
#'
#' @param u Numeric vector of log-abscissae.
#' @param y Numeric vector of log-probabilities.
#' @param inc Logical, TRUE for increasing (log-CDF), FALSE for decreasing (log-survival).
#' @return A list with strictly monotonic u and y.
#' @keywords internal
filter_strict_monotonic <- function(u, y, inc = FALSE) {
  ok <- which(is.finite(u) & is.finite(y) & (y > -300) & (y < 0))
  y_f <- y[ok]
  keep <- c(TRUE, if (inc) diff(cummax(y_f)) > 0 else diff(cummin(y_f)) < 0)
  list(u = u[ok][keep], y = y_f[keep])
}

#' Adaptive Log-Spline Builder with Richardson Error Control
#'
#' Constructs a Hyman monotonic cubic spline on u = log(x) with probability-space
#' interpolation error <= tol and zero wasted evaluations.
#'
#' @param u_min Numeric, lower log-bound.
#' @param u_max Numeric, upper log-bound.
#' @param eval_log_fn Function evaluating log probability on x = exp(u).
#' @param p_exp Numeric, initial grid power-law spacing exponent.
#' @param n_init Integer, initial knot count.
#' @param tol Numeric, probability-space error tolerance.
#' @param from_left Logical, whether power-law grid clusters to the left.
#' @param inc Logical, whether function is increasing.
#' @return A list with knots \code{u}, values \code{y}, and \code{n_evals}.
#' @keywords internal
build_adaptive_log_spline <- function(u_min, u_max, eval_log_fn, p_exp = 1.5,
                                      n_init = 51L, tol = 8e-13,
                                      from_left = FALSE, inc = FALSE) {
  s <- seq(0, 1, length.out = n_init + ((n_init + 1L) %% 2L))
  u <- u_min + (u_max - u_min) * (if (from_left) 1 - (1 - s)^p_exp else s^p_exp)
  y <- eval_log_fn(exp(u))
  
  for (pass in 1:6) {
    clean <- filter_strict_monotonic(u, y, inc = inc)
    u_ok <- clean$u
    y_ok <- clean$y
    n_ok <- length(u_ok)
    if (n_ok < 7L) break
    idx_even <- seq(1L, n_ok, by = 2L)
    if (idx_even[length(idx_even)] != n_ok) idx_even <- c(idx_even, n_ok)
    idx_odd <- seq(2L, n_ok - 1L, by = 2L)
    sp_even <- stats::splinefun(u_ok[idx_even], y_ok[idx_even], method = "hyman")
    y_pred_odd <- sp_even(u_ok[idx_odd])
    
    # Richardson Leave-Half-Out error estimate in probability space
    est_err <- abs(exp(y_pred_odd) - exp(y_ok[idx_odd])) * 0.0625 # / 16
    if (max(est_err, na.rm = TRUE) <= tol) break
    
    # 3-point max filter grades subdivision factors smoothly across neighbors
    n_odd <- length(est_err)
    err_pad <- c(est_err[1], est_err, est_err[n_odd])
    err_smooth <- pmax(est_err, err_pad[1:n_odd], err_pad[3:(n_odd + 2L)])
    sub_factors <- pmin((pass > 2L) * 4 + 4,
                        pmax(1, ceiling(sqrt(sqrt(err_smooth / tol)))))
    if (all(sub_factors <= 1)) break
    which_sub <- which(sub_factors > 1)
    new_u_list <- vector("list", length(which_sub))
    
    for (idx_k in seq_along(which_sub)) {
      j <- which_sub[idx_k]
      m_sub <- c(1L, sub_factors[j] + 1L)
      i0 <- idx_odd[j]
      u_b <- u_ok[i0]
      new_u_list[[idx_k]] <-
        c(seq(u_ok[i0 - 1L], u_b, length.out = m_sub[2])[-m_sub],
          seq(u_b, u_ok[min(n_ok, i0 + 1L)], length.out = m_sub[2])[-m_sub])
    }
    
    new_u <- sort(unique(unlist(new_u_list, use.names = FALSE)))
    
    # Proximity filter against existing knots
    if (length(new_u)) {
      fi <- findInterval(new_u, u)
      d_left <- new_u - u[pmax(1L, fi)]
      d_left[fi <= 0L] <- Inf
      d_right <- u[pmin(length(u), fi + 1L)] - new_u
      d_right[fi >= length(u)] <- Inf
      new_u <- new_u[pmin(d_left, d_right) > 1e-9]
    }
    
    if (!length(new_u)) break
    u <- c(u, new_u)
    y <- c(y, eval_log_fn(exp(new_u)))
    o <- order(u)
    u <- u[o]
    y <- y[o]
  }
  
  clean <- filter_strict_monotonic(u, y, inc = inc)
  list(u = clean$u, y = clean$y, n_evals = length(u))
}

#' Construct Adaptive Monotonic Spline Distribution Functions
#'
#' Factory function returning ultra-fast (\eqn{O(1)}) spline closures (\code{d}, \code{p},
#' \code{q}, \code{r}) for the convolution of arbitrary independent Gamma random variables.
#' Automatically selects between the 80-bit FFT Moschopoulos mixture and the saddlepoint
#' contour integral to evaluate knots, then fits adaptive Hyman monotonic cubic splines in
#' log-space with guaranteed error \eqn{\le} \code{tol}.
#'
#' @param shape Numeric vector of positive shape parameters.
#' @param rate Numeric vector of positive rate parameters.
#' @param scale Optional numeric vector of positive scale parameters (\code{1 / rate}).
#' @param nn Integer, quadrature node count for contour integration (default 40).
#' @param kmax Integer, crossover threshold for mixture vs. direct contour integration.
#'   Default automatically scales with component count: 2^21 for K <= 200, 2^22 for K <= 1200, 2^23 for K > 1200.
#' @param tol Numeric, maximum interpolation error tolerance in probability space (default 8e-13).
#' @param verbose Logical, whether to print engine selection information.
#' @return A list of closures of class \code{c("sumgamma_spline", "sumgamma_closure")}:
#' \describe{
#'   \item{\code{d(x, log = FALSE)}}{Density via analytic derivative of spline.}
#'   \item{\code{p(q, lower.tail = TRUE, log.p = FALSE)}}{Cumulative distribution / survival function.}
#'   \item{\code{q(p, lower.tail = TRUE, log.p = FALSE)}}{Quantiles via inverse spline evaluation.}
#'   \item{\code{r(n)}}{Random variate generator via inversion sampling.}
#' }
#' @export
sumgamma_spline <- function(shape, rate, scale = 1 / rate, nn = 40L,
                            kmax = NULL, tol = 8e-13, verbose = FALSE) {
  params <- prepare_gamma_params(shape, rate, scale)
  
  alpha <- params$shape
  tau <- params$rate
  b <- params$scale
  mu <- params$mu
  sigma <- params$sigma
  N_comp <- params$K
  sumalpha <- sum(alpha)
  
  if (is.null(kmax)) {
    kmax <- if (N_comp <= 200L) {
      2097152L
    } else if (N_comp <= 1200L) {
      4194304L
    } else {
      8388608L
    }
  }
  
  tau0 <- tau[N_comp]
  b0 <- b[1] # b is sorted descending, so b[1] is largest scale, b0 in code was b[N_comp]
  # In original code: b = 1/tau; b0 <- b[N_comp] (which was the scale corresponding to tau0 = tau[N_comp])
  # Since tau0 is tau[N_comp], its scale is 1/tau0
  scale0 <- 1 / tau0
  ab <- alpha * b
  a2 <- alpha * 0.5
  
  gq <- gauss_legendre_nodes(as.integer(nn))
  nodes_x <- gq$nodes
  nodes_w <- gq$weights
  
  psumgamma1 <- function(x, lower.tail = FALSE, log.p = TRUE) {
    psumgamma1_cpp(x, mu, as.integer(nn), alpha, b, ab, a2,
                   nodes_x, nodes_w, log.p, lower.tail)
  }
  
  # Test if required power-of-2 mixture size k is within kmax
  if (kmax >= 1 && psumgamma1((sumalpha + kmax) * scale0) < -35) {
    k <- as.numeric(kmax)
    while (k > 1 && psumgamma1((sumalpha + (k <- k * 0.5)) * scale0) < -35) {}
    k <- as.integer(2 * k)
    crossed <- FALSE
  } else {
    crossed <- TRUE
    k <- NA_integer_
  }
  
  if (verbose) {
    message(sprintf("[sumgamma_spline] N = %d | kmax = %g | Crossed: %s%s",
                    N_comp, kmax, crossed, if (!crossed) {
                      sprintf(" -> Using FFT Mixture (k = %d = 2^%d)", k, as.integer(log2(k)))
                    } else " -> Using Direct Contour Inversion"))
  }
  
  if (!crossed) {
    mix <- compute_mixture_weights_cpp(alpha, tau, k)
    w <- mix$w
    cw_lower <- mix$cw_lower
    cw_upper <- mix$cw_upper
    i0 <- mix$i0
    
    eval_log_lower <- function(x) {
      psumgamma_series_cpp(x, sumalpha, tau0, i0, w, cw_lower, cw_upper,
                           lower_tail = TRUE, log_p = TRUE)
    }
    eval_log_upper <- function(x) {
      psumgamma_series_cpp(x, sumalpha, tau0, i0, w, cw_lower, cw_upper,
                           lower_tail = FALSE, log_p = TRUE)
    }
  } else {
    eval_log_lower <- function(x) psumgamma1(x, lower.tail = TRUE, log.p = TRUE)
    eval_log_upper <- function(x) psumgamma1(x, lower.tail = FALSE, log.p = TRUE)
  }
  
  # Left-tail bracket
  sigma3_mu <- 3 * sigma / mu
  step_L <- exp(-max(0.02, min(0.4, sigma3_mu)))
  x_hi <- mu
  x_lo <- mu * step_L
  x_min <- 1e-250 * mu
  
  while (x_lo > x_min && eval_log_lower(x_lo) > -35) {
    x_hi <- x_lo
    x_lo <- x_lo * step_L
  }
  
  lo_u <- log(x_lo)
  hi_u <- log(x_hi)
  if (eval_log_lower(x_lo) <= -35) {
    for (. in 1:5) {
      mid_u <- 0.5 * (lo_u + hi_u)
      if (eval_log_lower(exp(mid_u)) > -35) hi_u <- mid_u else lo_u <- mid_u
    }
  }
  xL0 <- exp(0.5 * (lo_u + hi_u))
  
  # Right-tail bracket
  step_R <- exp(max(0.05, min(0.7, sigma3_mu)))
  x_lo <- mu
  x_hi <- mu * step_R
  xmax <- 1e250 * mu
  
  while (x_hi < xmax && eval_log_upper(x_hi) > -35) {
    x_lo <- x_hi
    x_hi <- x_hi * step_R
  }
  
  lo_u <- log(x_lo)
  hi_u <- log(x_hi)
  for (. in 1:5) {
    mid_u <- 0.5 * (lo_u + hi_u)
    if (eval_log_upper(exp(mid_u)) > -35) lo_u <- mid_u else hi_u <- mid_u
  }
  xR0 <- exp(0.5 * (lo_u + hi_u))
  
  # Power-law spacing exponents
  u_mu <- log(mu)
  u_L0 <- log(xL0)
  u_R0 <- log(xR0)
  p_R <- min(3.0, max(1.2, -log(log1p(sigma3_mu) / max(1e-4, u_R0 - u_mu)) * 0.721347520444482))
  p_L <- min(3.0, max(1.2, -log(log(1 + 3 * min(0.3 * mu, sigma) / mu) / max(1e-4, u_mu - u_L0)) * 0.721347520444482))
  
  fit_L <- build_adaptive_log_spline(u_L0, u_mu, eval_log_lower, p_exp = p_L,
                                     n_init = 51L, tol = tol, from_left = TRUE, inc = TRUE)
  fit_R <- build_adaptive_log_spline(u_mu, u_R0, eval_log_upper, p_exp = p_R,
                                     n_init = 51L, tol = tol, from_left = FALSE, inc = FALSE)
  
  # Share interior knots across u_mu to lock derivative
  nL <- length(fit_L$u) - 1L
  nL <- max(1L, nL - 1L):nL
  nR <- 2:min(3L, length(fit_R$u))
  clean0 <- filter_strict_monotonic(
    c(fit_L$u, fit_R$u[nR]),
    c(fit_L$y, log(-expm1(fit_R$y[nR]))),
    inc = TRUE
  )
  clean1 <- filter_strict_monotonic(
    c(fit_L$u[nL], fit_R$u),
    c(log(-expm1(fit_L$y[nL])), fit_R$y),
    inc = FALSE
  )
  
  p0 <- stats::splinefun(clean0$u, clean0$y, method = "hyman")
  p1 <- stats::splinefun(clean1$u, clean1$y, method = "hyman")
  neglogpmu <- -p0(u_mu)
  
  # Inverse quantile splines q0 and q1
  x00 <- seq(min(clean0$u), u_mu, length.out = 1001)
  x11 <- seq(u_mu, max(clean1$u), length.out = 1001)
  y00 <- -p0(x00)
  y11 <- -p1(x11)
  keep00 <- c(TRUE, diff(cummin(y00)) < 0)
  keep11 <- c(TRUE, diff(cummax(y11)) > 0)
  q0 <- stats::splinefun(y00[keep00], x00[keep00], method = "hyman")
  q1 <- stats::splinefun(y11[keep11], x11[keep11], method = "hyman")
  
  d_fn <- function(x, log = FALSE) {
    x_num <- as.numeric(x)
    i <- which(x_num < mu)
    if (log) {
      if (length(i)) {
        if (length(i) < length(x_num)) {
          d <- numeric(length(x_num))
          lx_i <- log(x_num[i]); lx_ni <- log(x_num[-i])
          d[i] <- log(p0(lx_i, 1) / x_num[i]) + p0(lx_i)
          d[-i] <- log(-p1(lx_ni, 1) / x_num[-i]) + p1(lx_ni)
          d
        } else {
          lx <- log(x_num)
          log(p0(lx, 1) / x_num) + p0(lx)
        }
      } else {
        lx <- log(x_num)
        log(-p1(lx, 1) / x_num) + p1(lx)
      }
    } else {
      if (length(i)) {
        if (length(i) < length(x_num)) {
          d <- numeric(length(x_num))
          lx_i <- log(x_num[i]); lx_ni <- log(x_num[-i])
          d[i] <- p0(lx_i, 1) * exp(p0(lx_i)) / x_num[i]
          d[-i] <- -p1(lx_ni, 1) * exp(p1(lx_ni)) / x_num[-i]
          d
        } else {
          lx <- log(x_num)
          p0(lx, 1) * exp(p0(lx)) / x_num
        }
      } else {
        lx <- log(x_num)
        -p1(lx, 1) * exp(p1(lx)) / x_num
      }
    }
  }
  
  p_fn <- function(q, lower.tail = TRUE, log.p = FALSE) {
    q_num <- as.numeric(q)
    i <- which(q_num < mu)
    if (lower.tail) {
      if (log.p) {
        if (length(i)) {
          if (length(i) < length(q_num)) {
            p <- numeric(length(q_num))
            p[i] <- p0(log(q_num[i]))
            p[-i] <- log(-expm1(p1(log(q_num[-i]))))
            p
          } else p0(log(q_num))
        } else log(-expm1(p1(log(q_num))))
      } else {
        if (length(i)) {
          if (length(i) < length(q_num)) {
            p <- numeric(length(q_num))
            p[i] <- exp(p0(log(q_num[i])))
            p[-i] <- -expm1(p1(log(q_num[-i])))
            p
          } else exp(p0(log(q_num)))
        } else -expm1(p1(log(q_num)))
      }
    } else {
      if (log.p) {
        if (length(i)) {
          if (length(i) < length(q_num)) {
            p <- numeric(length(q_num))
            p[i] <- log(-expm1(p0(log(q_num[i]))))
            p[-i] <- p1(log(q_num[-i]))
            p
          } else log(-expm1(p0(log(q_num))))
        } else p1(log(q_num))
      } else {
        if (length(i)) {
          if (length(i) < length(q_num)) {
            p <- numeric(length(q_num))
            p[i] <- -expm1(p0(log(q_num[i])))
            p[-i] <- exp(p1(log(q_num[-i])))
            p
          } else -expm1(p0(log(q_num)))
        } else exp(p1(log(q_num)))
      }
    }
  }
  
  q_fn <- function(p, lower.tail = TRUE, log.p = FALSE) {
    p_num <- as.numeric(p)
    if (lower.tail) {
      if (log.p) {
        i <- which(p_num < -neglogpmu)
        if (length(i)) {
          if (length(i) < length(p_num)) {
            q <- numeric(length(p_num))
            q[i] <- exp(q0(-p_num[i]))
            q[-i] <- exp(q1(-log(-expm1(p_num[-i]))))
            q
          } else exp(q0(-p_num))
        } else exp(q1(-log(-expm1(p_num))))
      } else {
        i <- which(p_num < exp(-neglogpmu))
        if (length(i)) {
          if (length(i) < length(p_num)) {
            q <- numeric(length(p_num))
            q[i] <- exp(q0(-log(p_num[i])))
            q[-i] <- exp(q1(-log1p(-p_num[-i])))
            q
          } else exp(q0(-log(p_num)))
        } else exp(q1(-log1p(-p_num)))
      }
    } else {
      if (log.p) {
        i <- which(p_num > log1p(-exp(-neglogpmu)))
        if (length(i)) {
          if (length(i) < length(p_num)) {
            q <- numeric(length(p_num))
            q[i] <- exp(q0(-log(-expm1(p_num[i]))))
            q[-i] <- exp(q1(-p_num[-i]))
            q
          } else exp(q0(-log(-expm1(p_num))))
        } else exp(q1(-p_num))
      } else {
        i <- which(p_num > -expm1(-neglogpmu))
        if (length(i)) {
          if (length(i) < length(p_num)) {
            q <- numeric(length(p_num))
            q[i] <- exp(q0(-log1p(-p_num[i])))
            q[-i] <- exp(q1(-log(p_num[-i])))
            q
          } else exp(q0(-log1p(-p_num)))
        } else exp(q1(-log(p_num)))
      }
    }
  }
  
  r_fn <- function(n) {
    r <- stats::rexp(as.integer(n))
    i <- which(r > neglogpmu)
    if (length(i)) {
      if (length(i) < n) {
        x <- numeric(n)
        x[i] <- exp(q0(r[i]))
        x[-i] <- exp(q1(-log1p(-exp(-r[-i]))))
        x
      } else exp(q0(r))
    } else exp(q1(-log1p(-exp(-r))))
  }
  
  out <- list(
    d = d_fn, p = p_fn, q = q_fn, r = r_fn,
    dsumgamma = d_fn, psumgamma = p_fn, qsumgamma = q_fn, rsumgamma = r_fn
  )
  attr(out, "params") <- params
  attr(out, "engine") <- "spline"
  attr(out, "k") <- k
  attr(out, "crossed") <- crossed
  attr(out, "bounds") <- c(xL0 = xL0, mu = mu, xR0 = xR0)
  attr(out, "n_knots") <- c(left = fit_L$n_evals, right = fit_R$n_evals)
  class(out) <- c("sumgamma_spline", "sumgamma_closure")
  out
}

