/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef KINETX_RUNGEKUTTA_H_
#define KINETX_RUNGEKUTTA_H_

#include "Kinetx/Network.h"

namespace Scine {
namespace Kinetx {
/**
 * @brief Base class for all Runge-Kutta methods/implementations
 */
class RungeKutta {
 public:
  /**
   * @brief Constructor
   * @param net The network of reactions.
   */
  RungeKutta(Network& net);

  /**
   * @brief Propagate the numerical integration by one time step.
   * @param concentrations The current concentrations. Updated inplace.
   * @param yFlux The current vertex (aggregate) flux. Updated inplace.
   * @param rFlux The current total edge (reaction) flux. Updated inplace.
   * @param rForwardFlux The current forward edge flux. Updated inplace.
   * @param rBackwardFlux The current backward edge flux. Updated inplace.
   * @param t The current time. Updated inplace.
   * @param dt The time increment. This may be updated inplace depending on the
   *           integration algorithm.
   */
  void propagate(Eigen::VectorXd& concentrations, Eigen::VectorXd& yFlux, Eigen::VectorXd& rFlux,
                 Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux, double& t, double& dt) const;
  /**
   * @brief Run the numerical integration.
   * @param y The input concentration.
   * @param tStart The start time.
   * @param dt The time increment.
   * @param rFlux The reaction edge flux (total).
   * @param rForwardFlux The forward reaction edge flux.
   * @param rBackwardFlux The backward reaction edge flux.
   * @param batchInterval The number of steps per batch.
   * @param nBatches The number of integration batches.
   * @param convergenceConcentrationChange The concentration convergence threshold.
   * @return The final concentrations, max. concentrations, and vertex fluxes as a matrix with columns
   *   in this order.
   */
  Eigen::MatrixXd runIntegration(Eigen::VectorXd y, double tStart, double dt, Eigen::VectorXd& rFlux,
                                 Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux,
                                 const unsigned int batchInterval = 1000, const unsigned int nBatches = 100000,
                                 const double convergenceConcentrationChange = 1e-10);

  Eigen::MatrixXd runIntegrationByTime(Eigen::VectorXd y, double tStart, double dt, Eigen::VectorXd& rFlux,
                                       Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux, const double tMax,
                                       const unsigned int batchInterval = 1000,
                                       const double convergenceConcentrationChange = 1e-10);

 protected:
  /**
   * @brief Keep track on the vertex and edge fluxes.
   * @param y The concentration after concentration propagation.
   * @param yInitial The original concentration.
   * @param yFlux The vertex flux. Updated inplace.
   * @param rFlux The total edge flux.
   * @param rForwardFlux The forward edge flux.
   * @param rBackwardFlux The backward edge flux.
   * @param dt The time increment.
   */
  void trackVertexAndEdgeFluxes(const Eigen::VectorXd& y, const Eigen::VectorXd& yInitial, Eigen::VectorXd& yFlux,
                                Eigen::VectorXd& rFlux, Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux,
                                const double& dt) const;
  /**
   * @brief Propagte the concentration.
   * @param concentrations The concentration.
   * @param t The current time.
   * @param dt The time increment.
   */
  virtual void propagateY(Eigen::VectorXd& concentrations, double& t, double& dt) const = 0;
  /**
   * @brief Calculate the concentration gradient.
   * @param concentrations The current concentration.
   * @return The concentration gradient.
   */
  Eigen::VectorXd g(const Eigen::VectorXd& concentrations) const;
  /**
   * @brief Calculate the concentration flux gradient.
   * @param concentrationsT0 The original concentraion.
   * @param concentrationsT1 The updated concentation.
   * @param flux The total reaction edge flux.
   * @param forwardFlux The reaction edge forward flux.
   * @param backwardFlux The reaction edge backward flux.
   * @return The vertex flux gradient.
   */
  Eigen::VectorXd gFlux(const Eigen::VectorXd& concentrationsT0, const Eigen::VectorXd& concentrationsT1,
                        Eigen::VectorXd& flux, Eigen::VectorXd& forwardFlux, Eigen::VectorXd& backwardFlux) const;
  /**
   * @brief Calculate the total edge flux gradient.
   * @param concentrationsT0 The original concentration.
   * @param concentrationsT1 The updated concentration.
   * @return The total edge flux.
   */
  Eigen::MatrixXd fFlux(const Eigen::VectorXd& concentrationsT0, const Eigen::VectorXd& concentrationsT1) const;
  /**
   * @brief Calculate the edge flux while keeping track on the direction.
   * @param concentrationsT0 The original concentration.
   * @param concentrationsT1 The updated concentration.
   * @return The total, forward, and backward flux (in this order) as a tuple.
   */
  std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> fFluxDirected(const Eigen::VectorXd& concentrationsT0,
                                                                              const Eigen::VectorXd& concentrationsT1) const;
  /**
   * Calculate the reaction rates f for each reaction using a mass-weight ansatz. This returns a matrix in case
   * several elementary steps are used for at least one reaction. Otherwise, this matrix will have only one column.
   * Hence, it will generally be densely populated.
   * @param concentrations The current concentrations.
   * @return The reaction rates as a vector/matrix.
   */
  Eigen::MatrixXd f(const Eigen::VectorXd& concentrations) const;
  /**
   * @brief The same as f(..) but keeping track on the directions.
   * @param concentrations The concentration.
   * @return The forward and backward reaction rates.
   */
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> fDirected(const Eigen::VectorXd& concentrations) const;
  Eigen::SparseMatrix<double> jacobi(const Eigen::VectorXd& concentrations) const;
  Network& _net;

  void printHeader(const Eigen::VectorXd& y, const double dt, const double t);
  bool printTimeAndCheckConvergenceStep(const Eigen::VectorXd& y, const Eigen::VectorXd& yOld, Eigen::VectorXd& yMax,
                                        const unsigned int batchInterval, const unsigned int iBatch, const double dt,
                                        const double t, const double convergenceConcentrationChange);
};

} /* namespace Kinetx */
} /* namespace Scine */

#endif // KINETX_RUNGEKUTTA_H_
