/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Includes */
#include <pybind11/pybind11.h>

void init_network(pybind11::module& m);
void init_network_builder(pybind11::module& m);
void init_reference_networks(pybind11::module& m);
void init_random_network_factory(pybind11::module& m);
void init_numerical_integration(pybind11::module& m);

PYBIND11_MODULE(scine_kinetx, m) {
  m.doc() = "Pybind11 Bindings for SCINE-Kinetx";

  // Ordering is important!
  init_network(m);
  init_network_builder(m);
  init_reference_networks(m);
  init_random_network_factory(m);
  init_numerical_integration(m);
}
