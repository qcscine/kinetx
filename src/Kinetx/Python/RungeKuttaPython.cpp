/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Include Internal Headers */
#include <Kinetx/RungeKutta/CashKarp5.h>     // numerical integration scheme.
#include <Kinetx/RungeKutta/ExplicitEuler.h> // numerical integration scheme.
#include <Kinetx/RungeKutta/ImplicitEuler.h> // numerical integration scheme.
#include <Kinetx/RungeKutta/RungeKutta.h>    // base class for the numerical integration.
/* Include Std and External Headers */
#include <pybind11/eigen.h>    // bind eigen3 objects.
#include <pybind11/pybind11.h> // python bindings.
#include <pybind11/stl.h>      // bind standard objects.
#include <iomanip>             // set_precision.
#include <memory>              // std::unique_ptr.

using namespace Scine::Kinetx;

/**
 * @brief Enum class for the different numerical integration schemes.
 */
enum class IntegratorName { CashKarp5, ExplicitEuler, ImplicitEuler };

/**
 * @brief A helper function to integrate the rate equations of a given reaction network in a batch-wise fashion.
 *        This function keeps track of the maximum concentration for each species in between integration batches.
 * @param network                         The reaction network.
 * @param yStart                          The starting concentrations.
 * @param tStart                          The starting time.
 * @param dt                              The time step (note that this may be ignored by some numerical
 *                                        integration schemes)
 * @param integratorName                  The enum for the integration scheme.
 * @param batchInterval                   The number of propagation steps per batch.
 * @param nBatches                        The maximum number of batches.
 * @param convergenceConcentrationChange  If the maximum change is concentration is lower than the given threshold. The
 *                                        integration is considered to be complete.
 * @param integrateByTime                 If true, the numerical integration is done up to the specified maxTime.
 * @param maxTime                         The maximum time to integrate. This is only used if integrateByTime is true.
 * @return A matrix of which the values of the first column are the final concentrations, the values of the second
 *         column are the maximum concentrations reached by the species, and the third column are the concentration
 *         flux through the species during the integration.
 */
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
integrateDifferentialEquations(Network& network, std::vector<double> yStart, const double tStart, const double dt,
                               const IntegratorName integratorName, const unsigned int batchInterval = 1000,
                               const unsigned int nBatches = 100000, const double convergenceConcentrationChange = 1e-10,
                               bool integrateByTime = false, double maxTime = 1e+5) {
  std::unique_ptr<RungeKutta> integrator;
  switch (integratorName) {
    case IntegratorName::CashKarp5:
      integrator = std::make_unique<CashKarp5>(network);
      break;
    case IntegratorName::ExplicitEuler:
      integrator = std::make_unique<ExplicitEuler>(network);
      break;
    case IntegratorName::ImplicitEuler:
      integrator = std::make_unique<ImplicitEuler>(network);
      break;
  }
  Eigen::VectorXd y = Eigen::Map<Eigen::VectorXd>(yStart.data(), yStart.size());
  Eigen::VectorXd rFlux, rBackwardsFlux, rForwardFlux;
  Eigen::MatrixXd concentrationResults;
  if (integrateByTime) {
    concentrationResults = integrator->runIntegrationByTime(y, tStart, dt, rFlux, rForwardFlux, rBackwardsFlux, maxTime,
                                                            batchInterval, convergenceConcentrationChange);
  }
  else {
    concentrationResults = integrator->runIntegration(y, tStart, dt, rFlux, rForwardFlux, rBackwardsFlux, batchInterval,
                                                      nBatches, convergenceConcentrationChange);
  }
  return std::make_tuple(concentrationResults, rFlux, rForwardFlux, rBackwardsFlux);
}

/**
 * @brief Resolve a string to its enum representation for the numerical integration schemes.
 * @param integrator The string.
 * @return The enum.
 */
IntegratorName resolveIntegratorName(std::string integrator) {
  if (integrator == "cash_karp_5") {
    return IntegratorName::CashKarp5;
  }
  else if (integrator == "explicit_euler") {
    return IntegratorName::ExplicitEuler;
  }
  else if (integrator == "implicit_euler") {
    return IntegratorName::ImplicitEuler;
  }
  throw std::runtime_error("The selected integrator is unknown.");
}

/**
 * @brief Python bindings for the numerical integration.
 * @param m
 */
void init_numerical_integration(pybind11::module& m) {
  pybind11::enum_<IntegratorName> integratorNames(m, "Integrator");
  integratorNames.value("cash_karp_5", IntegratorName::CashKarp5)
      .value("explicit_euler", IntegratorName::ExplicitEuler)
      .value("implicit_euler", IntegratorName::ImplicitEuler);

  m.def("get_integrator", &resolveIntegratorName, pybind11::arg("integrator"),
        R"delim(
      Resolve the integrator name string to its enum.
      :param integrator:   The name of the integrator to be used by kinetx.
    )delim");

  m.def("integrate", &integrateDifferentialEquations, pybind11::arg("network"), pybind11::arg("y_start"),
        pybind11::arg("t_start"), pybind11::arg("dt"), pybind11::arg("integrator_name"),
        pybind11::arg("batch_interval") = 1000, pybind11::arg("n_batches") = 100000,
        pybind11::arg("convergence") = 1e-10, pybind11::arg("integrateByTime") = false, pybind11::arg("maxTime") = 1e+5,
        R"delim(
      Propagate the time by n_batches x batch_interval x dt.
      :param network:         The reaction network.
      :param y_start:         The values to integrate (e.g. start concentrations).
      :param t_start:         The starting time.
      :param dt:              The time step.
      :param integrator_name: The integrator.
      :param batch_interval:  The number of time steps for each step batch (default: 1000).
      :param n_batches:       The number of batches (default: 100000).
      :param convergence:     The maximum change of concentration between two batches before considering
                              the kinetic model to be converged (default: 1e-10).
      :param integrateByTime  If true, the numerical integration is done up to the specified maxTime (default: False).
      :param maxTime          The maximum time to integrate. This is only used if integrateByTime is true (default: 1e+5).
    )delim");
}
