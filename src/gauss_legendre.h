#ifndef SUMGAMMA_GAUSS_LEGENDRE_H
#define SUMGAMMA_GAUSS_LEGENDRE_H

#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace sumgamma {

// Computes n-point Gauss-Legendre quadrature nodes and weights mapped to [0, 1].
// The nodes are strictly in (0, 1), sorted in ascending order.
// The weights sum to 2.0, adhering to the standard [-1, 1] weight convention
// so that: \int_0^1 f(x) dx = 0.5 * \sum_{i=1}^n w_i * f(x_i).
inline void compute_gauss_legendre_01(int n, std::vector<double>& nodes, std::vector<double>& weights) {
  nodes.resize(n);
  weights.resize(n);
  const int m = (n + 1) / 2;
  const double dn = static_cast<double>(n);
  
  for (int i = 0; i < m; ++i) {
    // Initial root approximation
    double z = std::cos(M_PI * (static_cast<double>(i) + 0.75) / (dn + 0.5));
    double pp = 0.0;
    
    // Newton-Raphson refinement on Legendre polynomial P_n(z)
    for (int iter = 0; iter < 100; ++iter) {
      double p1 = 1.0;
      double p2 = 0.0;
      for (int j = 1; j <= n; ++j) {
        double p3 = p2;
        p2 = p1;
        double dj = static_cast<double>(j);
        p1 = ((2.0 * dj - 1.0) * z * p2 - (dj - 1.0) * p3) / dj;
      }
      pp = dn * (z * p1 - p2) / (z * z - 1.0);
      double z1 = z;
      z = z1 - p1 / pp;
      if (std::abs(z - z1) < 1e-15) break;
    }
    
    double w = 2.0 / ((1.0 - z * z) * pp * pp);
    nodes[i] = 0.5 * (1.0 - z);
    nodes[n - 1 - i] = 0.5 * (1.0 + z);
    weights[i] = w;
    weights[n - 1 - i] = w;
  }
}

} // namespace sumgamma

#endif // SUMGAMMA_GAUSS_LEGENDRE_H

