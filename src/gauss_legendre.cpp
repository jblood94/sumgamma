#include <Rcpp.h>
#include "gauss_legendre.h"

using namespace Rcpp;

//' Gauss-Legendre Quadrature Nodes and Weights on [0, 1]
//'
//' Computes n-point Gauss-Legendre quadrature nodes and weights mapped to the unit interval [0, 1].
//'
//' @param n Integer, number of quadrature nodes (typically 20 to 60).
//' @return A list containing:
//' \describe{
//'   \item{nodes}{NumericVector of quadrature abscissae in (0, 1).}
//'   \item{weights}{NumericVector of weights (summing to 2.0, matching the [-1, 1] standard convention).}
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
List gauss_legendre_nodes(int n) {
  if (n < 1) stop("n must be at least 1");
  std::vector<double> nodes, weights;
  sumgamma::compute_gauss_legendre_01(n, nodes, weights);
  return List::create(
    Named("nodes") = wrap(nodes),
    Named("weights") = wrap(weights)
  );
}

