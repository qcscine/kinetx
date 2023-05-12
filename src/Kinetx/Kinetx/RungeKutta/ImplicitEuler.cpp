/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/RungeKutta/ImplicitEuler.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Scine {
namespace Kinetx {

ImplicitEuler::ImplicitEuler(Network& net) : RungeKutta(net) {
}

void ImplicitEuler::propagateY(Eigen::VectorXd& y, double& t, double& dt) const {
  auto j = this->jacobi(y);
  // i - (kroneckerProduct(a, i) * j
  // i - (kroneckerProduct(a, j)  ???
  Eigen::MatrixXd lhs = Eigen::MatrixXd::Identity(y.size(), y.size()) - dt * j;
  Eigen::VectorXd z = Eigen::VectorXd::Zero(y.size()); // y.size() * s
  Eigen::FullPivLU<Eigen::MatrixXd> lu(lhs);
  for (unsigned int newton = 0; newton < 100; newton++) {
    auto f = this->g(y + z); // do for s segments
    // i - (kroneckerProduct(a, f)  ???
    Eigen::VectorXd rhs = -z + dt * f;
    auto dz = lu.solve(rhs);
    z += dz;
  }
  t += dt;
  y += z;
}

} /* namespace Kinetx */
} /* namespace Scine */
