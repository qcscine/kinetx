__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import scine_kinetx as kx
import numpy as np


def test_numerical_integration_cash_karp_5():
    network_builder = kx.NetworkBuilder()
    n_compounds = 5
    n_reactions = 4
    n_channels_per_reaction = 1
    network_builder.reserve(n_compounds, n_reactions, n_channels_per_reaction)
    concentrations = [0.5, 0.4, 0.0, 0.0, 0.0]
    for i in range(n_compounds):
        network_builder.add_compound(1, str(i))
    network_builder.add_reaction([0.1], [0.05], [(0, 1)], [(2, 1)])
    network_builder.add_reaction([0.05], [0.05], [(1, 1)], [(3, 1)])
    network_builder.add_reaction(
        [0.02], [0.001], [(2, 1), (3, 1)], [(0, 1), (1, 1)])
    network_builder.add_reaction(
        [0.02], [0.0000001], [(0, 1), (1, 1)], [(4, 1)])
    network = network_builder.generate()
    solver = kx.Integrator.cash_karp_5
    concentration_data, r_flux, rf_flux, rb_flux = kx.integrate(network, np.asarray(
        concentrations), 0.0, 1e-7, solver, 1000, 5, 1e-10)
    reference_data = np.asarray(
        [3.337327e-02, 5.991047e-05, 6.674503e-02, 5.839152e-05, 3.998817e-01])
    reference_max = np.asarray([0.5, 0.4, 0.07103765, 0.00325947, 0.3998817])
    reference_flux = np.asarray(
        [1.85801505, 1.78337992, 1.45451083, 1.3798757, 0.40350422])
    difference = np.sum(np.abs(reference_data - concentration_data[:, 0]))

    assert difference < 1e-8
    difference = np.sum(np.abs(reference_max - concentration_data[:, 1]))
    assert difference < 1e-7
    difference = np.sum(np.abs(reference_flux - concentration_data[:, 2]))
    print(concentration_data[:, 2])
    assert difference < 1e-7

    concentration_data_2, r_flux_2, rf_flux_2, rb_flux_2 = kx.integrate(network, np.asarray(
        concentrations), 0.0, 1e-7, solver, 1000, 5, 1e-10, True, 1e+5)
    difference = np.sum(np.abs(concentration_data_2 - concentration_data))
    assert difference < 1e-9


def test_numerical_integration_implicit_euler():
    network_builder = kx.NetworkBuilder()
    n_compounds = 5
    n_reactions = 4
    n_channels_per_reaction = 1
    network_builder.reserve(n_compounds, n_reactions, n_channels_per_reaction)
    concentrations = [0.5, 0.4, 0.0, 0.0, 0.0]
    for i in range(n_compounds):
        network_builder.add_compound(1, str(i))
    network_builder.add_reaction([0.1], [0.05], [(0, 1)], [(2, 1)])
    network_builder.add_reaction([0.05], [0.05], [(1, 1)], [(3, 1)])
    network_builder.add_reaction(
        [0.02], [0.001], [(2, 1), (3, 1)], [(0, 1), (1, 1)])
    network_builder.add_reaction(
        [0.02], [0.0000001], [(0, 1), (1, 1)], [(4, 1)])
    network = network_builder.generate()
    solver = kx.Integrator.implicit_euler
    concentration_data, r_flux, rf_flux, rb_flux = kx.integrate(network, np.asarray(
        concentrations), 0.0, 7.75e+0, solver, 1000, 1000, 1e-10)
    reference_data = np.asarray(
        [3.337327e-02, 5.991047e-05, 6.674503e-02, 5.839152e-05, 3.998817e-01])
    reference_max = np.asarray([0.5, 0.4, 0.07103765, 0.00325947, 0.3998817])
    reference_flux = np.asarray(
        [1.85802001, 1.78338488, 1.45451295, 1.37987782, 0.40350706])
    difference = np.sum(np.abs(reference_data - concentration_data[:, 0]))
    assert difference < 1e-8
    difference = np.sum(np.abs(reference_max - concentration_data[:, 1]))
    assert difference < 1e-3
    difference = np.sum(np.abs(reference_flux - concentration_data[:, 2]))
    assert difference < 1


def test_numerical_integration_explicit_euler():
    network_builder = kx.NetworkBuilder()
    n_compounds = 5
    n_reactions = 4
    n_channels_per_reaction = 1
    network_builder.reserve(n_compounds, n_reactions, n_channels_per_reaction)
    concentrations = [0.5, 0.4, 0.0, 0.0, 0.0]
    for i in range(n_compounds):
        network_builder.add_compound(1, str(i))
    network_builder.add_reaction([0.1], [0.05], [(0, 1)], [(2, 1)])
    network_builder.add_reaction([0.05], [0.05], [(1, 1)], [(3, 1)])
    network_builder.add_reaction(
        [0.02], [0.001], [(2, 1), (3, 1)], [(0, 1), (1, 1)])
    network_builder.add_reaction(
        [0.02], [0.0000001], [(0, 1), (1, 1)], [(4, 1)])
    network = network_builder.generate()
    solver = kx.Integrator.explicit_euler
    concentration_data, r_flux, rf_flux, rb_flux = kx.integrate(network, np.asarray(
        concentrations), 0.0, 1e+0, solver, 5000, 10000, 1e-10)
    reference_data = np.asarray(
        [3.337327e-02, 5.991047e-05, 6.674503e-02, 5.839152e-05, 3.998817e-01])
    reference_max = np.asarray([0.5, 0.4, 0.07103765, 0.00325947, 0.3998817])
    reference_flux = np.asarray(
        [1.85802001, 1.78338488, 1.45451295, 1.37987782, 0.40350706])
    difference = np.sum(np.abs(reference_data - concentration_data[:, 0]))
    assert difference < 1e-7
    difference = np.sum(np.abs(reference_max - concentration_data[:, 1]))
    assert difference < 1
    difference = np.sum(np.abs(reference_flux - concentration_data[:, 2]))
    assert difference < 1
