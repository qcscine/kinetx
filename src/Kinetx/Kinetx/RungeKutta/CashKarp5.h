/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef KINETX_CASHKARP5_H_
#define KINETX_CASHKARP5_H_

#include "Kinetx/RungeKutta/RungeKutta.h"

namespace Scine {
namespace Kinetx {
/**
 * @brief Numerical integration according to the cash-karp-5 algorithm. This algorithm automatically
 *        generates the time step. (Ref. https://arxiv.org/pdf/1309.2710.pdf,
 *        https://doi.org/10.1016/j.jcp.2013.09.025)
 */
class CashKarp5 : public RungeKutta {
 public:
  /**
   * @brief Constructor
   * @param net The network of reactions.
   */
  CashKarp5(Network& net);

  /**
   * @brief
   */
  void propagateY(Eigen::VectorXd& y, double& t, double& dt) const override final;
};

} /* namespace Kinetx */
} /* namespace Scine */

#endif // KINETX_CASHKARP5_H_
