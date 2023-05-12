__copyright__ = """This file is part of SCINE Kinetx.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

from dev.conan.base import ScineConan


class ScineKinetxConan(ScineConan):
    name = "scine_kinetx"
    version = "2.0.0"
    url = "https://github.com/qcscine/kinetx"
    description = """ """
    options = {
        "shared": [True, False],
        "python": [True, False],
        "tests": [True, False],
        "coverage": [True, False],
        "microarch": ["detect", "none"]
    }
    default_options = {"shared": True, "python": False,
                       "tests": False, "coverage": False,
                       "microarch": "none"}
    exports = "dev/conan/*.py"
    exports_sources = ["dev/cmake/*", "src/*", "CMakeLists.txt", "README.rst",
                       "LICENSE.txt", "dev/conan/hook.cmake",
                       "dev/conan/glue/*"]
    requires = ["eigen/[~=3.3.7]"]
    cmake_name = "Kinetx"
    cmake_definitions = {
        "CMAKE_UNITY_BUILD": "ON",
        "CMAKE_UNITY_BUILD_BATCH_SIZE": 16
    }

    def package_info(self):
        super().package_info()
