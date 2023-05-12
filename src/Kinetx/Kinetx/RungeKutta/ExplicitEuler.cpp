/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/RungeKutta/ExplicitEuler.h"
#include <Eigen/Dense>

namespace Scine {
namespace Kinetx {

ExplicitEuler::ExplicitEuler(Network& net) : RungeKutta(net) {
}

void ExplicitEuler::propagateY(Eigen::VectorXd& y, double& t, double& dt) const {
  t += dt;
  y += dt * this->g(y);
}

} /* namespace Kinetx */
} /* namespace Scine */
