SCINE - KiNetX
==============

Introduction
------------

The primary objective of SCINE KiNetX is to model the kinetics of complex chemical reaction networks combined with rigorous uncertainty quantification. This is necessary to discover product distributions and reaction mechanisms and to evaluate their reliability. Besides routine kinetic analysis of arbitrary reaction networks, KiNetX can be utilized to steer the exploration of chemical reaction space to accelerate this exploration, by steering the exploration away from species which are kinetically irrelevant.

License and Copyright Information
---------------------------------

KiNetX is distributed under the BSD 3-clause "New" or "Revised" License.
For more license and copyright information, see the file ``LICENSE.txt`` in this
directory.

Installation and Usage
----------------------

The following software packages are required in order to compile SCINE KiNetX:

- A C++ compiler supporting the C++17 standard (GCC at least 7.3.0 or later)
- CMake (at least version 3.9.0)
- Eigen3 (at least version 3.3.2 or later)
- Boost (recommended: version 1.65.0 or later)

SCINE KiNetX can be built using a standard CMake/make setup::

    git submodule update --init
    mkdir build
    cd build
    cmake -DSCINE_BUILD_PYTHON_BINDINGS=ON ..
    make
    make test
    make install

How to Cite
-----------

When publishing results obtained with the SCINE database wrapper, please cite the corresponding
release as archived on Zenodo (please use the DOI of the respective release).

In addition, we kindly request you to cite the following articles when using KiNetX:

- \J. Proppe, M. Reiher, "Mechanism Deduction from Noisy Chemical Reaction Networks", *J. Chem. Theory Comput.*, **2019**, *15*, 357.
- \M. Bensberg, M. Reiher, "Concentration-Flux-Steered Mechanism Exploration with an Organocatalysis Application", *arXiv:2212.14135 [physics.chem-ph]*, **2022**.

Support and Contact
-------------------

In case you should encounter problems or bugs, please write a short message
to scine@phys.chem.ethz.ch.

Third-Party Libraries Used
--------------------------

SCINE KiNetX makes use of the following third-party libraries:

- `Boost <https://www.boost.org/>`_
- `Eigen <http://eigen.tuxfamily.org>`_
- `Google Test <https://github.com/google/googletest>`_
- `pybind11 <https://github.com/pybind/pybind11>`_
