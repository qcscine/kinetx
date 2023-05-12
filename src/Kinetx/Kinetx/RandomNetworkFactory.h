/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef KINETX_RANDOMNETWORKFACTORY_H_
#define KINETX_RANDOMNETWORKFACTORY_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Scine {
namespace Kinetx {

class Network;

/**
 * @brief A factory class for test networks.
 */
class RandomNetworkFactory {
 public:
  /**
   * @brief A random generator for physical networks
   *
   * The basis of this algorithm was called ``AutoNetGen`` in:
   *  "Mechanism Deduction from Noisy Chemical Reaction Networks"
   *   - Jonny Proppe and Markus Reiher
   *  https://doi.org/10.1021/acs.jctc.8b00310
   *
   * Small modifications may have been made w.r.t. to its original
   * implementation.
   *
   * @return Network A noisy test network.
   */
  static Network random();
};

} /* namespace Kinetx */
} /* namespace Scine */

#endif // KINETX_RANDOMNETWORKFACTORY_H_
