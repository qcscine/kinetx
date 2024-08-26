/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef KINETX_RUNGEKUTTA_H_
#define KINETX_RUNGEKUTTA_H_

#include "Kinetx/Integrator/Integrator.h"
#include "Kinetx/Network.h"

namespace Scine {
namespace Kinetx {
/**
 * @brief Base class for all Runge-Kutta methods/implementations
 */
class RungeKutta : public IntegratorBase {
 public:
  /**
   * @brief Constructor
   * @param net The network of reactions.
   */
  RungeKutta(Network& net);

  void propagate(Eigen::VectorXd& concentrations, Eigen::VectorXd& yFlux, Eigen::VectorXd& rFlux,
                 Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux, double& t, double& dt) const;

 protected:
  /**
   * @brief Propagate the concentration.
   * @param concentrations The concentration.
   * @param t The current time.
   * @param dt The time increment.
   */
  virtual void propagateY(Eigen::VectorXd& concentrations, double& t, double& dt) const = 0;
};

} /* namespace Kinetx */
} /* namespace Scine */

#endif // KINETX_RUNGEKUTTA_H_
