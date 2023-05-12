/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Include Internal Headers */
#include <Kinetx/Network.h>
/* Include Std and External Headers */
#include <pybind11/eigen.h>    // bind eigen3 objects.
#include <pybind11/pybind11.h> // python bindings.

using namespace Scine::Kinetx;

void init_network(pybind11::module& m) {
  pybind11::class_<Network> net(m, "Network");
  net.def(pybind11::init<Eigen::VectorXd, std::tuple<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>,
                         std::tuple<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int>>, std::vector<std::string>>(),
          R"delim(
        Constructor.
      :param masses:          The masses of the compounds in the network.
      :param rate_constants:  The rate constants ({forward, backward}) of all reactions and
                              channels (matrix order: (reaction index, channel index)).
      :param stoichiometry:   The stoichiometry of all reactions in the network ({lhs, rhs})
                              Matrix ordering is (reaction index, compound index).
      :param compound_labels: The compound labels.
    )delim");
  net.def_readonly("n_compounds", &Network::nCompounds,
                   R"delim(
      The number of compounds in the network.
    )delim");
  net.def_readonly("n_reactions", &Network::nReactions,
                   R"delim(
      The number of reactions in the network.
    )delim");
  net.def_readonly("masses", &Network::masses,
                   R"delim(
      The molecular masses of the compounds in the network.
    )delim");
  net.def_readonly("compound_labels", &Network::compoundLabels,
                   R"delim(
      Getter for the labels of all compounds in the network.
    )delim");
  net.def_readonly("rate_constants", &Network::rateConstants,
                   R"delim(
      The rate constants ({forward, backward}) of all reactions and
      channels (matrix order: (reaction index, channel index)).
    )delim");
  net.def_readonly("stoichiometry", &Network::stoichiometry,
                   R"delim(
      The stoichiometry of all reactions in the network ({lhs, rhs})
      Matrix ordering is (reaction index, compound index).
    )delim");
  net.def_readonly("total_stoichiometry", &Network::totalStoichiometry,
                   R"delim(
      The difference of backwards and forwards stoichiometry.
    )delim");
  net.def_readonly("added_stoichiometry", &Network::addedStoichiometry,
                   R"delim(
      The sum of forward and backwards stoichiometry.
    )delim");
}
