# sumgamma

[![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen.svg)]()
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**`sumgamma`** provides ultra-fast, high-precision distribution functions (CDF, PDF, quantile, and random variate generation) for the convolution of arbitrary independent Gamma random variables with non-identical shape and rate parameters:

$$X = \sum_{i=1}^K X_i, \quad X_i \sim \text{Gamma}(\text{shape} = \alpha_i, \text{rate} = \tau_i)$$

---

## Key Mathematical Architectures

1. **Saddlepoint Contour Integration (`sumgamma_integral`)**:
   - Shifts the Bromwich / Gil-Pelaez inversion contour through the exact saddlepoint $\hat{t}$ where $K'(\hat{t}) = x$.
   - Evaluates the lower tail ($x < \mu$) or upper tail ($x > \mu$) directly with zero numerical cancellation.
   - Singularity-adapted hyperbolic quadrature on $[0, z_0]$ resolving the nearest branch point.
   - Alternating tail cycles accelerated via **Wynn's $\varepsilon$-algorithm**.
   - $O(1)$ suffix-moment Taylor acceleration for distant scale components.

2. **80-Bit Extended-Precision FFT Mixture (`sumgamma_series`)**:
   - Fast Radix-2 FFT generator in 80-bit `long double` for the Moschopoulos (1985) infinite gamma mixture weights:
     $$f(x) = \sum_{k=0}^\infty w_k f_\Gamma(x; \rho + k, \tau_0)$$
   - Active Poisson-bracketed evaluation windows evaluated in native C++.
   - Exact two-stage categorical mixture random sampling (`$r(n)`).

3. **Adaptive Monotonic Spline Closures (`sumgamma_spline`)**:
   - Richardson leave-half-out adaptive Hyman monotonic cubic splines in $(\log x, \log p)$ space.
   - Guaranteed interpolation error $\le 10^{-12}$ with zero wasted evaluations.
   - Inverts the spline to construct $q(p)$, enabling sub-microsecond $O(1)$ quantile evaluation and inversion sampling.
   - Exact analytic derivative for density.

---

## Installation

You can install the development version of `sumgamma` from GitHub:

```r
# install.packages("remotes")
remotes::install_github("jblood94/sumgamma")
```

---

## Quickstart

### 1. Canonical Vectorized API

```r
library(sumgamma)

# Parameters for 3 independent Gamma components
shape <- c(1.5, 2.0, 3.2)
rate  <- c(0.8, 1.5, 3.0)

# Evaluate density, CDF, quantiles, and random variates
dsumgamma(c(2.5, 4.0), shape, rate)
psumgamma(c(2.5, 4.0), shape, rate)
qsumgamma(c(0.05, 0.5, 0.95), shape, rate)
rsumgamma(1000, shape, rate)
```

### 2. High-Performance Factory Closures

When evaluating the distribution repeatedly for fixed parameters, precomputing the distribution structure amortizes the setup cost:

```r
# Build fast O(1) spline closures (takes ~5 ms once)
dist_fn <- sumgamma_spline(shape, rate)

# Querying density, CDF, or quantiles takes nanoseconds
dist_fn$d(3.5)
dist_fn$p(3.5)
dist_fn$q(0.99)
dist_fn$r(100000)
```

### 3. S3 Distribution Interface

```r
dist <- sumgamma(shape, rate)

print(dist)
summary(dist)
plot(dist) # plots density and CDF side-by-side
```

### 4. Theoretical Moments

```r
sumgamma_moments(shape, rate)
# $mean: 4.273
# $var:  3.667
# $sd:   1.915
# $skewness: 0.956
# $kurtosis_excess: 1.482
```

---

## License

GPL (>= 3)
