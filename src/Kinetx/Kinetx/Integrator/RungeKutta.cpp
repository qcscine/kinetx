/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Include Class Header */
#include "Kinetx/Integrator/RungeKutta.h"
/* Include Std and External Headers */
#include <Eigen/Dense>  // Dense matrices
#include <Eigen/Sparse> // Sparse matrices

namespace Scine {
namespace Kinetx {

RungeKutta::RungeKutta(Network& net) : IntegratorBase(net) {
}

void RungeKutta::propagate(Eigen::VectorXd& concentrations, Eigen::VectorXd& yFlux, Eigen::VectorXd& rFlux,
                           Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux, double& t, double& dt) const {
  const Eigen::VectorXd yInitial = concentrations;
  this->propagateY(concentrations, t, dt);
  this->trackVertexAndEdgeFluxes(concentrations, yInitial, yFlux, rFlux, rForwardFlux, rBackwardFlux, dt);
}

} /* namespace Kinetx */
} /* namespace Scine */
