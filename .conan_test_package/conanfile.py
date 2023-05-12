__copyright__ = """This file is part of SCINE Kinetx.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

from conans import ConanFile, CMake


class TestPackageConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake_find_package"
    build_requires = "cmake/[>=3.18.0 <3.20.0 || >3.20.0]"
    exports_sources = "CMakeLists.txt", "test.cpp"

    def _configure(self):
        cmake = CMake(self)
        cmake.configure()
        return cmake

    def build(self):
        cmake = self._configure()
        cmake.build()

    def test(self):
        cmake = self._configure()
        cmake.test()

        if self.options["scine_kinetx"].python:
            self.output.info("Trying to import 'scine_kinetx'")
            import scine_kinetx
            self.output.info("Import worked!")
