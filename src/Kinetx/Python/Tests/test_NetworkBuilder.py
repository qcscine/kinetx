__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest

import scine_kinetx as kx


def test_constructor():
    nb = kx.NetworkBuilder()
    assert 1 == 1
