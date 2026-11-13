test_that("All three factory functions agree with single Gamma analytical case", {
  # When all rates are identical, sum of Gammas is exactly Gamma(sum(shape), rate)
  alpha_vec <- c(1.5, 2.5, 3.0)
  tau_common <- 2.0
  sum_a <- sum(alpha_vec) # 7.0
  
  fit_int <- sumgamma_integral(alpha_vec, rep(tau_common, 3))
  fit_ser <- sumgamma_series(alpha_vec, rep(tau_common, 3))
  fit_spl <- sumgamma_spline(alpha_vec, rep(tau_common, 3))
  
  x_vals <- c(1.0, 2.5, 3.5, 5.0, 7.0)
  true_p <- pgamma(x_vals, shape = sum_a, rate = tau_common)
  true_d <- dgamma(x_vals, shape = sum_a, rate = tau_common)
  
  # CDF checks
  expect_equal(fit_int$p(x_vals), true_p, tolerance = 1e-9)
  expect_equal(fit_ser$p(x_vals), true_p, tolerance = 1e-9)
  expect_equal(fit_spl$p(x_vals), true_p, tolerance = 1e-9)
  
  # PDF checks
  expect_equal(fit_int$d(x_vals), true_d, tolerance = 1e-9)
  expect_equal(fit_ser$d(x_vals), true_d, tolerance = 1e-9)
  expect_equal(fit_spl$d(x_vals), true_d, tolerance = 1e-9)
  
  # Quantile checks
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  true_q <- qgamma(probs, shape = sum_a, rate = tau_common)
  expect_equal(fit_int$q(probs), true_q, tolerance = 1e-8)
  expect_equal(fit_ser$q(probs), true_q, tolerance = 1e-8)
  expect_equal(fit_spl$q(probs), true_q, tolerance = 1e-8)
})

test_that("Spline builder achieves target tolerance <= 1e-11 on arbitrary non-identical rates", {
  set.seed(999)
  alpha <- c(0.8, 1.2, 2.5, 3.1)
  tau <- c(0.5, 1.5, 2.2, 4.0)
  
  fit_spl <- sumgamma_spline(alpha, tau, tol = 1e-12)
  fit_int <- sumgamma_integral(alpha, tau, nn = 40)
  
  mu <- sum(alpha / tau)
  x_test <- seq(0.2 * mu, 3.0 * mu, length.out = 25)
  
  # Compare spline CDF against saddlepoint contour integral ground truth
  err_cdf <- max(abs(fit_spl$p(x_test) - fit_int$p(x_test)))
  expect_lt(err_cdf, 1e-11)
  
  # Compare spline PDF against integral PDF
  err_pdf <- max(abs(fit_spl$d(x_test) - fit_int$d(x_test)))
  expect_lt(err_pdf, 1e-9)
})

