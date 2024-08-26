/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef KINETX_EXPLICITEULER_H_
#define KINETX_EXPLICITEULER_H_

#include "Kinetx/Integrator/RungeKutta.h"

namespace Scine {
namespace Kinetx {
/**
 * @brief Euler-type integration (dt * dy/dt).
 */
class ExplicitEuler : public RungeKutta {
 public:
  /**
   * @brief Constructor
   * @param net The network of reactions.
   */
  ExplicitEuler(Network& net);

  /**
   * @brief
   */
  void propagateY(Eigen::VectorXd& y, double& t, double& dt) const override final;
};

} /* namespace Kinetx */
} /* namespace Scine */

#endif // KINETX_EXPLICITEULER_H_
