/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/NetworkBuilder.h"
#include "Kinetx/Network.h"

namespace Scine {
namespace Kinetx {

Network NetworkBuilder::generate() {
  // Strip
  _ratesForward.conservativeResize(_nReactions, _nChannels);
  _ratesForward.prune(0.0);
  _ratesBackward.conservativeResize(_nReactions, _nChannels);
  _ratesBackward.prune(0.0);
  _stoichiometryForward.conservativeResize(_nReactions, _nCompounds);
  _stoichiometryForward.prune(0);
  _stoichiometryBackward.conservativeResize(_nReactions, _nCompounds);
  _stoichiometryBackward.prune(0);
  _masses.conservativeResize(_nCompounds);
  _labels.resize(_nCompounds);

  // Build Network
  Network net(_masses, std::tuple<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>{_ratesForward, _ratesBackward},
              std::tuple<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int>>{_stoichiometryForward, _stoichiometryBackward},
              _labels);
  return net;
}

void NetworkBuilder::reserve(unsigned int nCompounds, unsigned int nReactions, unsigned int nChannelsPerReaction) {
  if (nCompounds < _nCompounds) {
    throw std::runtime_error("Network Builder: Reserve requests less space than current content occupies.");
  }
  if (nReactions < _nReactions) {
    throw std::runtime_error("Network Builder: Reserve requests less space than current content occupies.");
  }
  if (nChannelsPerReaction < _nChannels) {
    throw std::runtime_error("Network Builder: Reserve requests less space than current content occupies.");
  }
  _ratesForward.conservativeResize(nReactions, nChannelsPerReaction);
  _ratesBackward.conservativeResize(nReactions, nChannelsPerReaction);
  _stoichiometryForward.conservativeResize(nReactions, nCompounds);
  _stoichiometryBackward.conservativeResize(nReactions, nCompounds);
  _masses.conservativeResize(nCompounds);
  _channelCounts.conservativeResize(nReactions);
  _channelCounts.segment(_nReactions, nReactions - _nReactions).setZero();
  _labels.resize(nCompounds);
}

void NetworkBuilder::addReaction(std::vector<double> rfs, std::vector<double> rbs,
                                 std::vector<std::pair<unsigned int, int>> lhs, std::vector<std::pair<unsigned int, int>> rhs) {
  // Error checks
  if (rfs.size() != rbs.size()) {
    throw std::runtime_error("Network Builder: The entered amount of forwards and backwards reaction rates differ.");
  }
  // Expand if required
  if (_nReactions + 1 > _ratesForward.rows()) {
    _ratesForward.conservativeResize(_nReactions + 1, std::max((int)rfs.size(), (int)_ratesForward.cols()));
    _ratesBackward.conservativeResize(_nReactions + 1, std::max((int)rfs.size(), (int)_ratesBackward.cols()));
    _stoichiometryForward.conservativeResize(_nReactions + 1, _stoichiometryForward.cols());
    _stoichiometryBackward.conservativeResize(_nReactions + 1, _stoichiometryBackward.cols());
    _channelCounts.conservativeResize(_nReactions + 1);
    _channelCounts(_nReactions) = 0;
  }
  if ((long int)rfs.size() > _ratesForward.cols()) {
    _ratesForward.conservativeResize(_ratesForward.rows(), rfs.size());
    _ratesBackward.conservativeResize(_ratesBackward.rows(), rfs.size());
  }
  // TODO: I replaced _nReactions + 1  by _nReactions. I do not understand why it was like this before...
  // Add data
  //  - Add Stoichiometry
  for (const auto& c : lhs) {
    if (c.first > _nCompounds - 1) {
      throw std::runtime_error("Network Builder: Reaction to be added references non existing compound.");
    }
    _stoichiometryForward.insert(_nReactions, c.first) = c.second;
  }
  for (const auto& c : rhs) {
    if (c.first > _nCompounds - 1) {
      throw std::runtime_error("Network Builder: Reaction to be added references non existing compound.");
    }
    _stoichiometryBackward.insert(_nReactions, c.first) = c.second;
  }
  //  - Add Rates
  for (unsigned int i = 0; i < rfs.size(); i++) {
    _ratesForward.insert(_nReactions, _channelCounts[_nReactions]) = rfs[i];
    _ratesBackward.insert(_nReactions, _channelCounts[_nReactions]) = rbs[i];
    _channelCounts[_nReactions] += 1;
  }
  _nChannels = std::max(_nChannels, (unsigned int)_channelCounts[_nReactions]);
  _nReactions++;
}

void NetworkBuilder::addReactionChannel(unsigned int reaction, double rf, double rb) {
  // Expand if required
  if (_channelCounts[reaction] + 1 > _ratesForward.cols()) {
    _ratesForward.conservativeResize(_ratesForward.rows(), _channelCounts[reaction] + 1);
    _ratesBackward.conservativeResize(_ratesBackward.rows(), _channelCounts[reaction] + 1);
  }
  // Add data
  _ratesForward.insert(reaction, _channelCounts[reaction]) = rf;
  _ratesBackward.insert(reaction, _channelCounts[reaction]) = rb;
  _channelCounts[reaction] += 1;
  _nChannels = std::max(_nChannels, (unsigned int)_channelCounts[reaction]);
}

void NetworkBuilder::addCompound(double mass, std::string label) {
  // Expand if required
  if (_nCompounds + 1 > _masses.size()) {
    _stoichiometryForward.conservativeResize(_stoichiometryForward.rows(), _nCompounds + 1);
    _stoichiometryBackward.conservativeResize(_stoichiometryBackward.rows(), _nCompounds + 1);
    _masses.conservativeResize(_nCompounds + 1);
    _labels.resize(_nCompounds + 1);
  }
  // Add data
  _masses[_nCompounds] = mass;
  _labels[_nCompounds] = label;
  _nCompounds++;
}

} /* namespace Kinetx */
} /* namespace Scine */
