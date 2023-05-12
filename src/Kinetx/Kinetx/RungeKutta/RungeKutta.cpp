/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Include Class Header */
#include "Kinetx/RungeKutta/RungeKutta.h"
/* Include Std and External Headers */
#include <Eigen/Dense>  // Dense matrices
#include <Eigen/Sparse> // Sparse matrices
#include <iomanip>      // set_precision.
#include <iostream>     // std::cout

namespace Scine {
namespace Kinetx {

RungeKutta::RungeKutta(Network& net) : _net(net) {
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> RungeKutta::fDirected(const Eigen::VectorXd& concentrations) const {
  /* Equation numbers given refer to the equations in:
   *  "Mechanism Deduction from Noisy Chemical Reaction Networks"
   *   - Jonny Proppe and Markus Reiher
   *  https://doi.org/10.1021/acs.jctc.8b00310
   */

  auto rf = std::get<0>(_net.rateConstants);
  auto rb = std::get<1>(_net.rateConstants);

  const auto sfT = std::get<0>(_net.stoichiometryTransposed);
  const auto sbT = std::get<1>(_net.stoichiometryTransposed);

  // Calculate forward and backward rates (Eq. 7)
  Eigen::VectorXd fRateBase = Eigen::VectorXd::Ones(_net.nReactions);
  Eigen::VectorXd bRateBase = Eigen::VectorXd::Ones(_net.nReactions);
  /* Parallelization is actually slower than serial execution for ca 140 compounds and reactions. */
  //#pragma omp parallel for
  for (unsigned int iRxn = 0; iRxn < _net.nReactions; ++iRxn) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(sfT, iRxn); it; ++it) {
      fRateBase.data()[iRxn] *= std::pow(concentrations[it.row()], it.value());
    }
    for (Eigen::SparseMatrix<int>::InnerIterator it(sbT, iRxn); it; ++it) {
      bRateBase.data()[iRxn] *= std::pow(concentrations[it.row()], it.value());
    }
  }
  rf.array().colwise() *= fRateBase.array();
  rb.array().colwise() *= bRateBase.array();
  return std::make_pair(rf, rb);
}

Eigen::MatrixXd RungeKutta::f(const Eigen::VectorXd& concentrations) const {
  const auto directedFluxes = fDirected(concentrations);
  // Calculate total rate (Eq. 8)
  return directedFluxes.first - directedFluxes.second;
}

Eigen::VectorXd RungeKutta::gFlux(const Eigen::VectorXd& concentrationsT0, const Eigen::VectorXd& concentrationsT1,
                                  Eigen::VectorXd& flux, Eigen::VectorXd& forwardFlux, Eigen::VectorXd& backwardFlux) const {
  const auto fDirectedT01 = fFluxDirected(concentrationsT0, concentrationsT1);
  const Eigen::MatrixXd cdt = _net.addedStoichiometry.cast<double>().transpose() * std::get<0>(fDirectedT01);
  // Row-wise sum down to a vector
  const Eigen::VectorXd vertexFlux = cdt.rowwise().sum();
  flux = std::get<0>(fDirectedT01).rowwise().sum();
  forwardFlux = std::get<1>(fDirectedT01).rowwise().sum();
  backwardFlux = std::get<2>(fDirectedT01).rowwise().sum();
  return vertexFlux;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
RungeKutta::fFluxDirected(const Eigen::VectorXd& concentrationsT0, const Eigen::VectorXd& concentrationsT1) const {
  const auto fT05Directed = fDirected(0.5 * (concentrationsT0 + concentrationsT1));
  const Eigen::MatrixXd& forwardFlux = fT05Directed.first;
  const Eigen::MatrixXd& backwardFlux = fT05Directed.second;
  const Eigen::MatrixXd totalFlux = (forwardFlux - backwardFlux).cwiseAbs();
  return std::make_tuple(totalFlux, forwardFlux, backwardFlux);
}

Eigen::MatrixXd RungeKutta::fFlux(const Eigen::VectorXd& concentrationsT0, const Eigen::VectorXd& concentrationsT1) const {
  // Eq. (19) from https://doi.org/10.1021/acs.jctc.8b00310
  const auto fT05 = f(0.5 * (concentrationsT0 + concentrationsT1));
  // Approximate f(0.5 (t_u + t_u-1)) as f(0.5 (t_u + t_u-1)) ~ 0.5 (f(t_u) + f(t_u-1))
  const Eigen::MatrixXd fT01 = fT05.cwiseAbs();
  return fT01;
}

Eigen::VectorXd RungeKutta::g(const Eigen::VectorXd& concentrations) const {
  /* Equation numbers given refer to the equations in:
   *  "Mechanism Deduction from Noisy Chemical Reaction Networks"
   *   - Jonny Proppe and Markus Reiher
   *  https://doi.org/10.1021/acs.jctc.8b00310
   */
  // Calculate the total reaction rates (Eq. 8).
  const auto totalRates = f(concentrations);
  // Calculate time dependent concentrations (Eq. 9)
  const Eigen::MatrixXd cdt = _net.totalStoichiometry.cast<double>().transpose() * totalRates;
  const Eigen::VectorXd g = cdt.rowwise().sum();
  return g;
}

Eigen::SparseMatrix<double> RungeKutta::jacobi(const Eigen::VectorXd& concentrations) const {
  const auto& sf = std::get<0>(_net.stoichiometry);
  const auto& sb = std::get<1>(_net.stoichiometry);
  // Contract rate constants per channel into totals
  const auto& rfs = std::get<0>(_net.rateConstants);
  const auto& rbs = std::get<1>(_net.rateConstants);
  const Eigen::VectorXd rf = rfs.rowwise().sum(); // * Eigen::VectorXd::Ones(rfs.cols()); // verbose rowwise sum
  const Eigen::VectorXd rb = rbs.rowwise().sum(); // * Eigen::VectorXd::Ones(rbs.cols()); // verbose rowwise sum

  // Calculate derivative
  Eigen::SparseMatrix<double> fRateDeriv(_net.nReactions, _net.nCompounds);
  Eigen::SparseMatrix<double> bRateDeriv(_net.nReactions, _net.nCompounds);
  // TODO OMP? + Use transposed stoichiometry matrices for more efficient traversing through the matrices
  std::vector<Eigen::Triplet<double>> fRateTriplets;
  fRateTriplets.reserve(sf.nonZeros());
  for (int col = 0; col < sf.outerSize(); ++col) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(sf, col); it; ++it) {
      fRateTriplets.push_back(Eigen::Triplet<double>(it.row(), col, 1.0));
    }
  }
  fRateDeriv.setFromTriplets(fRateTriplets.begin(), fRateTriplets.end());
  for (int col = 0; col < sf.outerSize(); ++col) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(sf, col); it; ++it) {
      double tmp = fRateDeriv.coeff(it.row(), col);
      fRateDeriv.col(col) *= std::pow(concentrations[col], it.value());
      fRateDeriv.coeffRef(it.row(), col) = tmp * it.value() * std::pow(concentrations[col], it.value() - 1);
    }
  }

  std::vector<Eigen::Triplet<double>> bRateTriplets;
  bRateDeriv.reserve(sb.nonZeros());
  for (int col = 0; col < sb.outerSize(); ++col) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(sb, col); it; ++it) {
      bRateTriplets.push_back(Eigen::Triplet<double>(it.row(), col, 1.0));
    }
  }
  bRateDeriv.setFromTriplets(bRateTriplets.begin(), bRateTriplets.end());
  for (int col = 0; col < sb.outerSize(); ++col) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(sb, col); it; ++it) {
      double tmp = bRateDeriv.coeff(it.row(), col);
      bRateDeriv.col(col) *= std::pow(concentrations[col], it.value());
      bRateDeriv.coeffRef(it.row(), col) = tmp * it.value() * std::pow(concentrations[col], it.value() - 1);
    }
  }
  for (unsigned int col = 0; col < _net.nCompounds; ++col) {
    fRateDeriv.col(col) = fRateDeriv.col(col).cwiseProduct(rf);
    bRateDeriv.col(col) = bRateDeriv.col(col).cwiseProduct(rb);
  }
  return (_net.totalStoichiometry.cast<double>().transpose() * (fRateDeriv - bRateDeriv)).transpose();
}

Eigen::MatrixXd RungeKutta::runIntegration(Eigen::VectorXd y, double t, double dt, Eigen::VectorXd& rFlux,
                                           Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux,
                                           const unsigned int batchInterval, const unsigned int nBatches,
                                           const double convergenceConcentrationChange) {
  Eigen::VectorXd yMax = y;
  Eigen::VectorXd yInt = Eigen::VectorXd::Zero(y.rows());
  rFlux = Eigen::VectorXd::Zero(this->_net.nReactions);
  rForwardFlux = Eigen::VectorXd::Zero(this->_net.nReactions);
  rBackwardFlux = Eigen::VectorXd::Zero(this->_net.nReactions);
  printHeader(y, dt, t);
  for (unsigned int iBatch = 0; iBatch < nBatches; ++iBatch) {
    const Eigen::VectorXd yOld = y;
    // Propagate and keep track of the concentration flux.
    for (unsigned int i = 0; i < batchInterval; ++i) {
      this->propagate(y, yInt, rFlux, rForwardFlux, rBackwardFlux, t, dt);
    } // for i
    if (printTimeAndCheckConvergenceStep(y, yOld, yMax, batchInterval, iBatch, dt, t, convergenceConcentrationChange))
      break;
  } // for iBatch
  Eigen::MatrixXd toReturn(yMax.rows(), 3);
  toReturn.col(0) = y;
  toReturn.col(1) = yMax;
  toReturn.col(2) = yInt;
  return toReturn;
}

void RungeKutta::printHeader(const Eigen::VectorXd& y, const double dt, const double t) {
  std::cout << "#             t/s          "
            << " dt/s         "
            << " concentrations   mass_x_con_check_sum  max-change" << std::endl;
  std::cout << "0             " << std::scientific << std::setprecision(6) << t << "  " << dt << "  " << y.transpose()
            << "  " << y.dot(this->_net.masses) << std::endl;
}

bool RungeKutta::printTimeAndCheckConvergenceStep(const Eigen::VectorXd& y, const Eigen::VectorXd& yOld,
                                                  Eigen::VectorXd& yMax, const unsigned int batchInterval,
                                                  const unsigned int iBatch, const double dt, const double t,
                                                  const double convergenceConcentrationChange) {
  yMax.array() = yMax.array().max(y.array());
  const double maxChange = (y - yOld).array().abs().maxCoeff();
  std::cout << std::scientific << std::setprecision(6) << (double)(iBatch + 1) * batchInterval << "  " << t << "  "
            << dt << "  " << y.transpose() << "  " << y.dot(this->_net.masses) << "  " << maxChange << std::endl;
  if (maxChange != maxChange)
    throw std::runtime_error("NaN detected during numerical integration.");
  if (maxChange < convergenceConcentrationChange) {
    std::cout << "converged better than " << convergenceConcentrationChange << std::endl;
    return true;
  }
  return false;
}

Eigen::MatrixXd RungeKutta::runIntegrationByTime(Eigen::VectorXd y, double t, double dt, Eigen::VectorXd& rFlux,
                                                 Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux,
                                                 const double tMax, const unsigned int batchInterval,
                                                 const double convergenceConcentrationChange) {
  Eigen::VectorXd yMax = y;
  Eigen::VectorXd yInt = Eigen::VectorXd::Zero(y.rows());
  rFlux = Eigen::VectorXd::Zero(this->_net.nReactions);
  rForwardFlux = Eigen::VectorXd::Zero(this->_net.nReactions);
  rBackwardFlux = Eigen::VectorXd::Zero(this->_net.nReactions);
  printHeader(y, dt, t);
  unsigned int iBatch = 0;
  while (t < tMax) {
    const Eigen::VectorXd yOld = y;
    // Propagate and keep track of the concentration flux.
    for (unsigned int i = 0; i < batchInterval; ++i) {
      this->propagate(y, yInt, rFlux, rForwardFlux, rBackwardFlux, t, dt);
    } // for i
    if (printTimeAndCheckConvergenceStep(y, yOld, yMax, batchInterval, iBatch, dt, t, convergenceConcentrationChange))
      break;
    iBatch++;
  } // while t < tMax
  Eigen::MatrixXd toReturn(yMax.rows(), 3);
  toReturn.col(0) = y;
  toReturn.col(1) = yMax;
  toReturn.col(2) = yInt;
  return toReturn;
}

void RungeKutta::propagate(Eigen::VectorXd& concentrations, Eigen::VectorXd& yFlux, Eigen::VectorXd& rFlux,
                           Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux, double& t, double& dt) const {
  const Eigen::VectorXd yInitial = concentrations;
  this->propagateY(concentrations, t, dt);
  this->trackVertexAndEdgeFluxes(concentrations, yInitial, yFlux, rFlux, rForwardFlux, rBackwardFlux, dt);
}

void RungeKutta::trackVertexAndEdgeFluxes(const Eigen::VectorXd& y, const Eigen::VectorXd& yInitial,
                                          Eigen::VectorXd& yFlux, Eigen::VectorXd& rFlux, Eigen::VectorXd& rForwardFlux,
                                          Eigen::VectorXd& rBackwardFlux, const double& dt) const {
  Eigen::VectorXd flux, forwardFlux, backwardFlux;
  yFlux += dt * this->gFlux(yInitial, y, flux, forwardFlux, backwardFlux);
  rFlux += dt * flux;
  rForwardFlux += dt * forwardFlux;
  rBackwardFlux += dt * backwardFlux;
}

} /* namespace Kinetx */
} /* namespace Scine */
