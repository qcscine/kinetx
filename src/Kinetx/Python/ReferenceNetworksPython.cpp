/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kinetx/Network.h>
#include <Kinetx/ReferenceNetworks.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Kinetx::ReferenceNetworks;

void init_reference_networks(pybind11::module& m) {
  pybind11::module ref = m.def_submodule("reference_networks");
  ref.doc() = R"(The ``reference_networks`` submodule defines a set of known fixed networks.)";
  ref.def("get_bray_liebhafsky", &getBrayLiebhafsky,
          R"delim(
      Get a model for the Bray–Liebhafsky reaction.
      The model used here can be found in:
      "Improvement of the stoichiometric network analysis for determination of
      instability conditions of complex nonlinear reaction systems"
      Ljiljana Kolar-Anic, Zeljko Cupic, Guy Schmitz, Slobodan Anic
      Chemical Engineering Science, 65, (2010), 3718-3728
      https://www.sciencedirect.com/science/article/pii/S0009250910001569

      :return The newtork and a set of initial concentrations.
    )delim");
}
