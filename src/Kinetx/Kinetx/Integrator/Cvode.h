/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef KINETX_CVODE_H_
#define KINETX_CVODE_H_

#include "Kinetx/Integrator/Integrator.h"
#include "Kinetx/Network.h"
#include <memory>

namespace Scine {
namespace Kinetx {
/**
 * @brief Sundials implementation called CVODE
 */
class Cvode : public Integrator {
 public:
  /**
   * @brief Constructor
   * @param net The network of reactions.
   */
  Cvode(Network& net);

  /// @bried Destructor.
  ~Cvode();

  Eigen::MatrixXd runIntegrationByTime(Eigen::VectorXd y, double t, double dt, Eigen::VectorXd& rFlux,
                                       Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux, const double tMax,
                                       const unsigned int batchInterval, const double convergenceConcentrationChange);

  Eigen::MatrixXd runIntegration(Eigen::VectorXd y, double t, double dt, Eigen::VectorXd& rFlux, Eigen::VectorXd& rForwardFlux,
                                 Eigen::VectorXd& rBackwardFlux, const unsigned int batchInterval,
                                 const unsigned int nBatches, const double convergenceConcentrationChange);

  void propagate(Eigen::VectorXd& concentrations, Eigen::VectorXd& yFlux, Eigen::VectorXd& rFlux,
                 Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux, double& t, double& dt) const;

 private:
  class Impl;
  std::unique_ptr<Impl> _pimpl;
};

} /* namespace Kinetx */
} /* namespace Scine */

#endif // KINETX_CVODE_H_