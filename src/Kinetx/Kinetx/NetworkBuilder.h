/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef KINETX_NETWORKBUILDER_H_
#define KINETX_NETWORKBUILDER_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Scine {
namespace Kinetx {

class Network;

/**
 * @brief A class allowing for easier building of reaction networks.
 */
class NetworkBuilder {
 public:
  /// @brief Default constructor.
  NetworkBuilder() = default;
  /**
   * @brief Generates the final network.
   * @return Network The network containing all previously added data.
   */
  Network generate();
  /**
   * @brief Reserves space in the underlying data objects, allowing for faster
   *        inserts.
   * @param nCompounds           The maximum number of compounds.
   * @param nReactions           The maximum number of reactions.
   * @param nChannelsPerReaction The maximum number of reaction channels per reaction.
   */
  void reserve(unsigned int nCompounds, unsigned int nReactions, unsigned int nChannelsPerReaction);
  /**
   * @brief Adds a new reaction to the network, auto expands fields if required.
   * @param rfs Reaction rates for the forward reaction (one per channel).
   * @param rbs Reaction rates for the backward reaction (one per channel).
   * @param lhs Stoichiometry of the LHS of the reaction.
   *            (Format: {{Compound1, Equivalents1}, {Compound1, Equivalents2}, ...}).
   * @param rhs Stoichiometry of the RHS of the reaction.
   *            (Format: {{Compound1, Equivalents1}, {Compound1, Equivalents2}, ...}).
   */
  void addReaction(std::vector<double> rfs, std::vector<double> rbs, std::vector<std::pair<unsigned int, int>> lhs,
                   std::vector<std::pair<unsigned int, int>> rhs);
  /**
   * @brief Adds a single reaction channel to an existing reaction.
   *
   * (Expands fields if required.)
   *
   * @param reaction The number (index, 0 based) of the reaction to add to.
   * @param rf       The reaction rate for the forward reaction.
   * @param rb       The reaction rate fot the backward reaction.
   */
  void addReactionChannel(unsigned int reaction, double rf, double rb);
  /**
   * @brief Adds a single new compound to the network.
   *
   * (Expands fields if required.)
   *
   * @param mass  The molecular mass of the compound to add.
   * @param label Optional: The label of the compound.
   */
  void addCompound(double mass, std::string label = "");

 private:
  unsigned int _nCompounds = 0;
  unsigned int _nReactions = 0;
  unsigned int _nChannels = 0;
  Eigen::VectorXi _channelCounts;
  Eigen::VectorXd _masses;
  std::vector<std::string> _labels;
  Eigen::SparseMatrix<double> _ratesForward;
  Eigen::SparseMatrix<double> _ratesBackward;
  Eigen::SparseMatrix<int> _stoichiometryForward;
  Eigen::SparseMatrix<int> _stoichiometryBackward;
};

} /* namespace Kinetx */
} /* namespace Scine */

#endif // KINETX_NETWORKBUILDER_H_
