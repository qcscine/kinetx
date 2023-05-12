/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kinetx/Network.h>
#include <Kinetx/NetworkBuilder.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Kinetx;

void init_network_builder(pybind11::module& m) {
  pybind11::class_<NetworkBuilder> builder(m, "NetworkBuilder");

  builder.def(pybind11::init<>(),
              R"delim(
      Initialize a network builder.
    )delim");
  builder.def("generate", &NetworkBuilder::generate,
              R"delim(
      Generates the final network.
    )delim");
  builder.def("reserve", &NetworkBuilder::reserve, pybind11::arg("n_compounds"), pybind11::arg("n_reactions"),
              pybind11::arg("n_channels_per_reaction"),
              R"delim(
      Reserves space in the underlying data objects, allowing for faster inserts.
      :param n_compounds: The maximum number of compounds.
      :param n_reactions: The maximum number of reactions.
      :param n_channels_per_reaction: The maximum number of reaction channels per reaction.
    )delim");
  builder.def("add_reaction", &NetworkBuilder::addReaction, pybind11::arg("forward_rates"),
              pybind11::arg("backward_rates"), pybind11::arg("left_hand_side"), pybind11::arg("right_hand_side"),
              R"delim(
      Adds a new reaction to the network, auto expands fields if required.
      :param forward_rates: Reaction rates for the forward reaction (one per channel).
      :param backward_rates: Reaction rates for the backward reaction (one per channel).
      :param left_hand_side: Stoichiometry of the LHS of the reaction.
                             (Format: [(Compound1, Equivalents1), (Compound2, Equivalents2), ...] ).
      :param right_hand_side: Stoichiometry of the RHS of the reaction.
                             (Format: [(Compound1, Equivalents1), (Compound2, Equivalents2), ...] ).
    )delim");
  builder.def("add_reaction_channel", &NetworkBuilder::addReactionChannel, pybind11::arg("reaction"),
              pybind11::arg("forward_rate"), pybind11::arg("backward_rate"),
              R"delim(
      Adds a single reaction channel to an existing reaction.
      (Expands fields if required.)
      :param reaction: The number (index, 0 based) of the reaction to add to.
      :param forward_rate: The reaction rate for the forward reaction.
      :param backward_rate: The reaction rate fot the backward reaction.
    )delim");
  builder.def("add_compound", &NetworkBuilder::addCompound, pybind11::arg("mass"), pybind11::arg("label") = "",
              R"delim(
      Adds a single new compound to the network.
      (Expands fields if required.)
      :param mass: The molecular mass of the compound to add.
      :param label: Optional: The label of the compound.
    )delim");
}
