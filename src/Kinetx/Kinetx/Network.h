/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef KINETX_NETWORK_H_
#define KINETX_NETWORK_H_

#include "Kinetx/Network.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Scine {
namespace Kinetx {

class Network {
 public:
  /**
   * @brief Constructor
   *
   * @param masses        The masses of all compounds
   * @param rateConstants The reaction rates for all reactions.
   *                      Data layout:
   *                      Rates: {forward, backward}
   *                      Matrix order: (reaction index, channel index)
   * @param stoichiometry The stoichiometry tuple for all reactions.
   *                      Data layout:
   *                      Stoichiometry: {lhs, rhs, total}
   *                      Matrix order: (reaction index, compound index)
   * @param labels        The label of all compounds in the network.
   */
  Network(Eigen::VectorXd masses, std::tuple<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> rateConstants,
          std::tuple<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int>> stoichiometry, std::vector<std::string> labels = {});

  // void load(std::string path);
  // void write(std::string path);

  /// @brief The number of compounds in the network.
  const unsigned int nCompounds;
  /// @brief The number of reactions in the network.
  const unsigned int nReactions;
  /// @brief The molecular masses of the compounds in the network.
  const Eigen::VectorXd masses;
  /**
   * @brief The rate constants ({forward, backward}) of all reactions and
   *        channels (matrix order: (reaction index, channel index)).
   */
  const std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> rateConstants;
  /**
   * @brief The stoichiometry of all reactions in the network ({lhs, rhs, lhs - rhs, lhs + rhs})
   *
   * The total stoichiometry is given in forward direction (rhs - lhs).
   * Matrix ordering is (reaction index, compound index).
   */
  const std::tuple<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int>> stoichiometry;

  const std::tuple<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int>> stoichiometryTransposed;
  ///@brief The difference of backwards and forwards stoichiometry.
  const Eigen::SparseMatrix<int> totalStoichiometry;
  ///@brief The sum of forward and backwards stoichiometry.
  const Eigen::SparseMatrix<int> addedStoichiometry;
  /// @brief The label of all compounds in the network.
  const std::vector<std::string> compoundLabels;

 private:
};

} /* namespace Kinetx */
} /* namespace Scine */

#endif // KINETX_NETWORK_H_
