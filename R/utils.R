#' Validate and Prepare Gamma Convolution Parameters
#'
#' Validates shape and rate/scale parameter vectors, aggregates identical rates,
#' and returns sorted parameters.
#'
#' @param shape Numeric vector of positive shape parameters (alpha).
#' @param rate Numeric vector of positive rate parameters (tau).
#' @param scale Optional numeric vector of positive scale parameters (1 / rate).
#' @return A list containing:
#' \describe{
#'   \item{shape}{Numeric vector of aggregated shape parameters sorted by rate ascending.}
#'   \item{rate}{Numeric vector of distinct rate parameters sorted ascending.}
#'   \item{scale}{Numeric vector of scale parameters (1 / rate) sorted descending.}
#'   \item{mu}{Theoretical mean of the sum.}
#'   \item{var}{Theoretical variance of the sum.}
#'   \item{sigma}{Theoretical standard deviation of the sum.}
#'   \item{K}{Number of distinct gamma components after aggregation.}
#' }
#' @keywords internal
#' @export
prepare_gamma_params <- function(shape, rate, scale = 1 / rate) {
  if (missing(rate) && !missing(scale)) {
    rate <- 1 / scale
  } else if (!missing(rate) && missing(scale)) {
    scale <- 1 / rate
  }
  
  if (!is.numeric(shape) || !is.numeric(rate)) {
    stop("'shape' and 'rate' must be numeric vectors")
  }
  
  if (length(shape) != length(rate)) {
    if (length(shape) == 1L) {
      shape <- rep.int(shape, length(rate))
    } else if (length(rate) == 1L) {
      rate <- rep.int(rate, length(shape))
    } else {
      stop("Lengths of 'shape' and 'rate' must match")
    }
  }
  
  if (any(shape <= 0 | !is.finite(shape))) {
    stop("All 'shape' parameters must be finite and strictly positive")
  }
  if (any(rate <= 0 | !is.finite(rate))) {
    stop("All 'rate' parameters must be finite and strictly positive")
  }
  
  # Aggregate components that share identical rate parameters
  if (anyDuplicated(rate)) {
    agg <- rowsum(shape, rate, reorder = TRUE)
    rate <- as.numeric(rownames(agg))
    shape <- as.numeric(agg)
  } else {
    o <- order(rate)
    rate <- rate[o]
    shape <- shape[o]
  }
  
  scale <- 1 / rate
  mu <- sum(shape * scale)
  var_val <- sum(shape * scale^2)
  sigma <- sqrt(var_val)
  
  list(
    shape = shape,
    rate = rate,
    scale = scale,
    mu = mu,
    var = var_val,
    sigma = sigma,
    K = length(shape)
  )
}

