
/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef KINETX_REFERENCENETWORK_H_
#define KINETX_REFERENCENETWORK_H_

#include <Eigen/Dense>
#include <utility>

namespace Scine {
namespace Kinetx {
class Network;
namespace ReferenceNetworks {
/**
 * @brief Get a model for the Bray–Liebhafsky reaction.
 *
 * The model used here can be found in:
 * "Improvement of the stoichiometric network analysis for determination of
 * instability conditions of complex nonlinear reaction systems"
 * Ljiljana Kolar-Anic, Zeljko Cupic, Guy Schmitz, Slobodan Anic;
 * Chemical Engineering Science, 65, (2010), 3718-3728
 * https://www.sciencedirect.com/science/article/pii/S0009250910001569
 *
 * @return The newtork and a set of initial concentrations.
 */
std::pair<Network, Eigen::VectorXd> getBrayLiebhafsky();

} /* namespace ReferenceNetworks */
} /* namespace Kinetx */
} /* namespace Scine */

#endif // KINETX_REFERENCENETWORK_H_
