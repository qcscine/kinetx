/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/ReferenceNetworks.h"
#include "Kinetx/Network.h"

namespace Scine {
namespace Kinetx {
namespace ReferenceNetworks {

std::pair<Network, Eigen::VectorXd> getBrayLiebhafsky() {
  // Bray–Liebhafsky
  // Ref: https://www.sciencedirect.com/science/article/pii/S0009250910001569
  Eigen::VectorXd masses(6);
  masses << 126.90, 143.91, 159.91, 269.81, 253.81, 34.01;
  double j0 = 0.028;
  Eigen::SparseMatrix<double> rf(15, 1);
  rf.insert(0, 0) = 1.375e2;
  rf.insert(1, 0) = 4.79e10;
  rf.insert(2, 0) = 5.00e3;
  rf.insert(3, 0) = 3.0e11;
  rf.insert(4, 0) = 1.486e4;
  rf.insert(5, 0) = 5.00e5;
  rf.insert(7, 0) = 4.4606e-5;
  rf.insert(8, 0) = 0.155 * j0;
  rf.insert(9, 0) = j0;
  rf.insert(10, 0) = j0;
  rf.insert(11, 0) = j0;
  rf.insert(12, 0) = j0;
  rf.insert(13, 0) = j0;
  rf.insert(14, 0) = j0;
  rf /= 60;
  Eigen::SparseMatrix<double> rb(15, 1);
  rb.insert(0, 0) = 7.91e7;
  rb.insert(2, 0) = 3.15e8;
  rb.insert(3, 0) = 46.97;
  rb /= 60;
  Eigen::SparseMatrix<int> sf(15, 6);
  sf.insert(0, 0) = 1;
  sf.insert(1, 0) = 1;
  sf.insert(1, 2) = 1;
  sf.insert(2, 3) = 1;
  sf.insert(3, 0) = 1;
  sf.insert(3, 1) = 1;
  sf.insert(4, 1) = 1;
  sf.insert(4, 5) = 1;
  sf.insert(5, 3) = 1;
  sf.insert(5, 5) = 1;
  sf.insert(6, 2) = 1;
  sf.insert(6, 5) = 1;
  sf.insert(7, 5) = 1;
  sf.insert(9, 4) = 1;
  sf.insert(10, 5) = 1;
  sf.insert(11, 0) = 1;
  sf.insert(12, 1) = 1;
  sf.insert(13, 2) = 1;
  sf.insert(14, 3) = 1;
  Eigen::SparseMatrix<int> sb(15, 6);
  sb.insert(0, 1) = 1;
  sb.insert(0, 2) = 1;
  sb.insert(1, 3) = 1;
  sb.insert(2, 1) = 2;
  sb.insert(3, 4) = 1;
  sb.insert(4, 0) = 1;
  sb.insert(5, 1) = 1;
  sb.insert(5, 2) = 1;
  sb.insert(7, 2) = 1;
  sb.insert(8, 5) = 1;
  Network net(masses, std::tuple<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>{rf, rb},
              std::tuple<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int>>{sf, sb});
  Eigen::VectorXd concentrations = Eigen::VectorXd::Zero(net.nCompounds);
  concentrations[0] = 1.70e-8;
  concentrations[1] = 9.20e-8;
  concentrations[2] = 3.20e-7;
  concentrations[3] = 5.30e-10;
  concentrations[4] = 1.00e-5;
  concentrations[5] = 1.55e-1;
  return {net, concentrations};
}

} /* namespace ReferenceNetworks */
} /* namespace Kinetx */
} /* namespace Scine */
