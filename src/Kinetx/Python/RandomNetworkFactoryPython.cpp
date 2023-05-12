/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kinetx/Network.h>
#include <Kinetx/RandomNetworkFactory.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Kinetx;

void init_random_network_factory(pybind11::module& m) {
  pybind11::class_<RandomNetworkFactory> factory(m, "RandomNetworkFactory");

  factory.def("random", &RandomNetworkFactory::random,
              R"delim(
      A random generator for physical networks
      The basis of this algorithm was called ``AutoNetGen`` in:
       "Mechanism Deduction from Noisy Chemical Reaction Networks"
        - Jonny Proppe and Markus Reiher
       https://doi.org/10.1021/acs.jctc.8b00310

      Small modifications may have been made w.r.t. to its original
      implementation.

      :return Network A noisy test network.
    )delim");
}
