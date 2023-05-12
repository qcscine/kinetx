/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/Network.h"
#include "Kinetx/RandomNetworkFactory.h"
#include "Kinetx/ReferenceNetworks.h"
#include "Kinetx/RungeKutta/CashKarp5.h"
#include "Kinetx/RungeKutta/ExplicitEuler.h"
#include "Kinetx/RungeKutta/ImplicitEuler.h"
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace Scine::Kinetx;

int main() {
  auto ret = ReferenceNetworks::getBrayLiebhafsky();
  auto net = ret.first;
  auto concentrations = ret.second;

  // Minimal Bromate Oscillator
  // https://pubs.acs.org/doi/pdf/10.1021/j100259a030
  //  Eigen::VectorXd masses(10);
  //  masses << 127.90, 79.90, 1.01, 112.91, 96.91, 159.81, 18.02, 111.90, 140.12, 140.12;
  //  Eigen::SparseMatrix<double> rf(7, 1);
  //  rf.insert(0, 0) = 2.1;
  //  rf.insert(1, 0) = 2.0e9;
  //  rf.insert(2, 0) = 8.0e9;
  //  rf.insert(3, 0) = 1.0e4;
  //  rf.insert(4, 0) = 6.5e5;
  //  rf.insert(5, 0) = 9.6;
  //  rf.insert(6, 0) = 4.0e7;
  //  Eigen::SparseMatrix<double> rb(7, 1);
  //  rb.insert(0, 0) = 1.0e4;
  //  rb.insert(1, 0) = 5.0e-5;
  //  rb.insert(2, 0) = 110.0;
  //  rb.insert(3, 0) = 2.0e7;
  //  rb.insert(4, 0) = 2.4e7;
  //  rb.insert(5, 0) = 1.3e-4;
  //  rb.insert(6, 0) = 2.1e-10;
  //  Eigen::SparseMatrix<int> sf(7, 10);
  //  sf.insert(0, 0) = 1;
  //  sf.insert(0, 1) = 1;
  //  sf.insert(0, 2) = 2;
  //  sf.insert(1, 1) = 1;
  //  sf.insert(1, 2) = 1;
  //  sf.insert(1, 3) = 1;
  //  sf.insert(2, 1) = 1;
  //  sf.insert(2, 2) = 1;
  //  sf.insert(2, 4) = 1;
  //  sf.insert(3, 0) = 1;
  //  sf.insert(3, 2) = 1;
  //  sf.insert(3, 3) = 1;
  //  sf.insert(4, 2) = 1;
  //  sf.insert(4, 7) = 1;
  //  sf.insert(4, 8) = 1;
  //  sf.insert(5, 6) = 1;
  //  sf.insert(5, 7) = 1;
  //  sf.insert(5, 9) = 1;
  //  sf.insert(6, 3) = 2;
  //  Eigen::SparseMatrix<int> sb(7, 10);
  //  sb.insert(0, 3) = 1;
  //  sb.insert(0, 4) = 1;
  //  sb.insert(1, 4) = 2;
  //  sb.insert(2, 5) = 1;
  //  sb.insert(2, 6) = 1;
  //  sb.insert(3, 6) = 1;
  //  sb.insert(3, 7) = 2;
  //  sb.insert(4, 3) = 1;
  //  sb.insert(4, 9) = 1;
  //  sb.insert(5, 8) = 1;
  //  sb.insert(5, 0) = 1;
  //  sb.insert(5, 2) = 2;
  //  sb.insert(6, 0) = 1;
  //  sb.insert(6, 2) = 1;
  //  sb.insert(6, 4) = 1;
  //  Network net(masses, {rf, rb}, {sf, sb, sb - sf});
  //  Eigen::VectorXd concentrations = Eigen::VectorXd::Zero(net.nCompounds);
  //  concentrations[0] = 6.0e-2;
  //  concentrations[1] = 3.0e-4;
  //  concentrations[2] = 1.5;
  //  concentrations[6] = 22.0;
  //  concentrations[8] = 1.5e-4;
  // concentrations[0] = 0.1999e-1;
  // concentrations[1] = 0.3347e-6;
  // concentrations[2] = 0.1500e+1;
  // concentrations[3] = 0.4606e-10;
  // concentrations[4] = 0.8271e-6;
  // concentrations[5] = 0.3020e-4;
  // concentrations[6] = 0.0;
  // concentrations[7] = 0.2627e-9;
  // concentrations[8] = 0.1427e-3;
  // concentrations[9] = 0.7321e-5;

  // Random Network
  // auto net = NetworkFactory::random();
  // Eigen::VectorXd concentrations(net.nCompounds);
  // concentrations.setZero();
  // concentrations[0] = 0.3;

  // Explicit Euler
  // ExplicitEuler solver(net);
  // Cash-Karp
  CashKarp5 solver(net);
  // Implicit Euler
  // ImplicitEuler solver(net);
  double t = 0.0;
  double dt = 1e-8;
  // for (unsigned int i=0;i<1e5;i++){
  std::cout << std::scientific << std::setprecision(6) << t << "  " << dt << "  " << concentrations.transpose() << "  "
            << concentrations.dot(net.masses) << std::endl;
  Eigen::VectorXd concentrationFlux = Eigen::VectorXd::Zero(concentrations.size());
  Eigen::VectorXd edgeFlux = Eigen::VectorXd::Zero(net.nReactions);
  Eigen::VectorXd forwardEdgeFlux = Eigen::VectorXd::Zero(net.nReactions);
  Eigen::VectorXd backwardEdgeFlux = Eigen::VectorXd::Zero(net.nReactions);
  for (unsigned int i = 0; i < 2e8; i++) {
    double told = t;
    solver.propagate(concentrations, concentrationFlux, edgeFlux, forwardEdgeFlux, backwardEdgeFlux, t, dt);
    if (t != told && i % 1000 == 0)
      std::cout << std::scientific << std::setprecision(6) << t << "  " << dt << "  " << concentrations.transpose()
                << "  " << concentrations.dot(net.masses) << std::endl;
  }

  return 1;
}
