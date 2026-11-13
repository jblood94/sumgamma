test_that("Canonical d/p/q/r functions and S3 object work as expected", {
  alpha <- c(1.2, 2.3)
  tau <- c(1.0, 3.0)
  
  # Canonical functions
  d_vals <- dsumgamma(c(1, 2), alpha, tau)
  p_vals <- psumgamma(c(1, 2), alpha, tau)
  q_vals <- qsumgamma(c(0.1, 0.9), alpha, tau)
  r_vals <- rsumgamma(100, alpha, tau)
  
  expect_equal(length(d_vals), 2)
  expect_equal(length(p_vals), 2)
  expect_equal(length(q_vals), 2)
  expect_equal(length(r_vals), 100)
  expect_true(all(d_vals > 0))
  expect_true(all(p_vals >= 0 & p_vals <= 1))
  expect_true(all(q_vals > 0))
  expect_true(all(r_vals > 0))
  
  # Direct simulation
  r_dir <- rsumgamma_direct(500, alpha, tau)
  expect_equal(length(r_dir), 500)
  expect_true(all(r_dir > 0))
  
  # Moments
  moms <- sumgamma_moments(alpha, tau)
  expect_equal(moms$mean, sum(alpha / tau))
  expect_equal(moms$var, sum(alpha / tau^2))
  
  # Characteristic function: phi(0) == 1
  expect_equal(cf_sumgamma(0, alpha, tau), complex(real = 1, imaginary = 0))
  
  # Diagnostics
  diag <- sumgamma_diagnose(alpha, tau)
  expect_equal(diag$components, 2)
  expect_type(diag$recommended_engine, "character")
  
  # S3 Object
  dist <- sumgamma(alpha, tau)
  expect_s3_class(dist, "sumgamma")
  expect_output(print(dist), "Sum of 2 Independent Gamma")
  
  s <- summary(dist)
  expect_s3_class(s, "summary.sumgamma")
  expect_output(print(s), "Summary of Sum of 2 Gamma")
})

test_that("Parameter validation catches invalid inputs", {
  expect_error(prepare_gamma_params(c(1, 2, 3), c(1, 2)))
  expect_error(prepare_gamma_params(c(-1, 2), c(1, 2)))
  expect_error(prepare_gamma_params(c(1, 2), c(0, 2)))
  expect_error(prepare_gamma_params("a", 1))
})
