/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/Integrator/CashKarp5.h"
#include "Kinetx/Integrator/Cvode.h"
#include "Kinetx/Integrator/ExplicitEuler.h"
#include "Kinetx/Integrator/ImplicitEuler.h"
#include "Kinetx/Network.h"
#include "Kinetx/RandomNetworkFactory.h"
#include "Kinetx/ReferenceNetworks.h"
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace Scine::Kinetx;

int main() {
  auto ret = ReferenceNetworks::getBrayLiebhafsky();
  auto net = ret.first;
  auto concentrations = ret.second;

  // Explicit Euler
  // ExplicitEuler solver(net);
  // Implicit Euler
  // ImplicitEuler solver(net);
  // Cash-Karp
  // CashKarp5 solver(net);
  // Cvode
  Cvode solver(net);
  double t = 0.0;
  double dt = 1e-1;
  Eigen::VectorXd edgeFlux = Eigen::VectorXd::Zero(net.nReactions);
  Eigen::VectorXd forwardEdgeFlux = Eigen::VectorXd::Zero(net.nReactions);
  Eigen::VectorXd backwardEdgeFlux = Eigen::VectorXd::Zero(net.nReactions);
  for (unsigned int i = 0; i < 2e6; i++) {
    double told = t;
    Eigen::VectorXd concentrationFlux = Eigen::VectorXd::Zero(concentrations.size());
    solver.propagate(concentrations, concentrationFlux, edgeFlux, forwardEdgeFlux, backwardEdgeFlux, t, dt);
    if (t != told && i % 1000 == 0)
      std::cout << std::scientific << std::setprecision(6) << t << "  " << dt << " | " << concentrations.transpose()
                << " | " << concentrations.dot(net.masses) << " | " << (concentrationFlux / dt).transpose() << std::endl;
  }

  return 1;
}
