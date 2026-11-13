test_that("Gauss-Legendre nodes match statmod::gauss.quad to machine precision", {
  skip_if_not_installed("statmod")
  for (n in c(10, 20, 40, 60)) {
    gq_ref <- statmod::gauss.quad(n)
    ref_nodes <- (gq_ref$nodes + 1) * 0.5
    ref_weights <- gq_ref$weights
    
    res <- gauss_legendre_nodes(n)
    expect_equal(sort(res$nodes), sort(ref_nodes), tolerance = 1e-14)
    expect_equal(sort(res$weights), sort(ref_weights), tolerance = 1e-14)
  }
})

test_that("Moschopoulos series and Saddlepoint Contour Integral agree to high precision", {
  set.seed(123)
  N <- 4
  tau <- sort(rexp(N, 1/N), decreasing = FALSE) # sorted ascending
  alpha <- rexp(N, N) * tau
  sumalpha <- sum(alpha)
  b <- 1 / tau
  ab <- alpha * b
  a2 <- alpha * 0.5
  mu <- sum(ab)
  
  # Compute mixture weights
  mix <- compute_mixture_weights_cpp(alpha, tau, M = 2048)
  
  # Test points across left tail, center, and right tail
  x_eval <- c(0.1 * mu, 0.5 * mu, mu, 1.8 * mu, 3.5 * mu)
  
  # CDF comparison
  p_series <- psumgamma_series_cpp(x_eval, sumalpha, mix$tau0, mix$i0,
                                   mix$w, mix$cw_lower, mix$cw_upper,
                                   lower_tail = TRUE, log_p = FALSE)
  
  # Saddlepoint requires scale b sorted descending:
  o_desc <- order(b, decreasing = TRUE)
  p_integral <- psumgamma1_cpp(x_eval, mu, 40,
                               alpha[o_desc], b[o_desc], ab[o_desc], a2[o_desc],
                               numeric(0), numeric(0),
                               log_p = FALSE, lower_tail = TRUE)
  
  expect_equal(p_series, p_integral, tolerance = 1e-10)
  
  # PDF comparison
  d_series <- dsumgamma_series_cpp(x_eval, sumalpha, mix$tau0, mix$i0,
                                   mix$w, log_d = FALSE)
  d_integral <- dsumgamma1_cpp(x_eval, mu, 40,
                               alpha[o_desc], b[o_desc], ab[o_desc], a2[o_desc],
                               numeric(0), numeric(0),
                               log_d = FALSE)
  
  expect_equal(d_series, d_integral, tolerance = 1e-10)
})

test_that("rsumgamma_series_cpp produces samples with matching theoretical mean and variance", {
  set.seed(456)
  alpha <- c(2.0, 3.5, 1.2)
  tau <- c(1.5, 2.0, 4.0)
  sumalpha <- sum(alpha)
  b <- 1 / tau
  mu_theo <- sum(alpha * b)
  var_theo <- sum(alpha * b^2)
  
  mix <- compute_mixture_weights_cpp(alpha, tau, M = 1024)
  samples <- rsumgamma_series_cpp(50000, sumalpha, mix$tau0, mix$i0, mix$cw_lower)
  
  expect_equal(mean(samples), mu_theo, tolerance = 0.02)
  expect_equal(var(samples), var_theo, tolerance = 0.05)
})

