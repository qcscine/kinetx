/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef KINETX_IMPLICITEULER_H_
#define KINETX_IMPLICITEULER_H_

#include "Kinetx/RungeKutta/RungeKutta.h"

namespace Scine {
namespace Kinetx {
/**
 * @brief Implicit Euler-type integration algorithm (y_k+1 = dt f(t_k+1, y_k+1)),
 *        where f(t_k+1, y_k+1) is the reaction rate and k the integration step.
 *        Requires the Jacobian of the reaction system.
 */
class ImplicitEuler : public RungeKutta {
 public:
  /**
   * @brief Constructor
   * @param net The network of reactions.
   */
  ImplicitEuler(Network& net);

  /**
   * @brief
   */
  void propagateY(Eigen::VectorXd& y, double& t, double& dt) const override final;
};

} /* namespace Kinetx */
} /* namespace Scine */

#endif // KINETX_IMPLICITEULER_H_
