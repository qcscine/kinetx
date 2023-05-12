# Load source file lists
include(Files.cmake)

if(SCINE_PARALLELIZE)
  find_package(OpenMP)
endif()

# Kinetx library
add_library(Kinetx ${KINETX_LIB_FILES})

target_include_directories(Kinetx PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/Scine>
)
set_target_properties(Kinetx PROPERTIES OUTPUT_NAME kinetx)
if(NOT BUILD_SHARED_LIBS)
  set_target_properties(Kinetx PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
target_compile_options(Kinetx PRIVATE
  $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>
)
target_compile_definitions(Kinetx
 PRIVATE
  $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
)
include(AddEigen)
add_eigen(Kinetx PUBLIC)
target_link_libraries(Kinetx
  PRIVATE
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

# Add namespaced aliases
add_library(Scine::Kinetx ALIAS Kinetx)

# -- Define Install
# Headers
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/
  DESTINATION include/Scine/Kinetx/
  FILES_MATCHING PATTERN "*.h"
)
# Targets
install(TARGETS Kinetx
  EXPORT kinetxTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Component-level cmake package files
scine_install_component_cmake_files(
  COMPONENT  ScineKinetx
  EXPORT_NAME kinetxTargets
)

# Tests
if(SCINE_BUILD_TESTS)
  include(ImportGTest)
  import_gtest()
  # Gtest
  add_executable(Kinetx_tests ${KINETX_TEST_FILES})
  target_link_libraries(Kinetx_tests PRIVATE
    GTest::Main
    GMock::GMock
    Scine::Kinetx
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
    ${CMAKE_DL_LIBS}
  )
  add_test(NAME Kinetx COMMAND Kinetx_tests)

  # # App tests
  # find_package(PythonInterp REQUIRED)
  # add_test(
  #   NAME Sparrow_App
  #   COMMAND ${PYTHON_EXECUTABLE} -B -m pytest test_sparrow.py${TEST_SELECTION} --junitxml=${CMAKE_CURRENT_BINARY_DIR}/app_report.xml
  #   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Tests/AppTests
  # )
  # set_tests_properties(Sparrow_App PROPERTIES
  #   ENVIRONMENT SCINE_MODULE_PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{SCINE_MODULE_PATH}
  #   ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PYTHONPATH}
  #   ENVIRONMENT PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PATH}
  # )
endif()

# Set the RPATH to be used when installing.
if(APPLE)
  set(CMAKE_INSTALL_RPATH "@loader_path/../lib")
else()
  # $ORIGIN refers to the installation directory of the binary.
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()

# Executable
add_executable(KinetxApp ${KINETX_APP_FILES})
add_executable(Scine::KinetxApp ALIAS KinetxApp)

set_target_properties(KinetxApp PROPERTIES OUTPUT_NAME kinetx)
target_link_libraries(KinetxApp PRIVATE
  Scine::Kinetx
  ${CMAKE_DL_LIBS}
)
install(TARGETS KinetxApp RUNTIME DESTINATION bin)

# Add link library specific to MinGW
if(WIN32 AND MINGW)
  target_link_libraries(KinetxApp PRIVATE ws2_32)
endif()


# Python Bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  include(ImportPybind11)
  import_pybind11()
  set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})

  # Python module
  pybind11_add_module(scine_kinetx ${KINETX_PYTHON_CPPS})
  set_target_properties(scine_kinetx PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/scine_kinetx
  )
  target_include_directories(scine_kinetx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Python)
  target_link_libraries(scine_kinetx PRIVATE Scine::Kinetx)
  install(CODE
    "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade ${CMAKE_CURRENT_BINARY_DIR}
                     RESULT_VARIABLE retcode)
     if(NOT \$\{retcode\} EQUAL 0)
       message(FATAL_ERROR \"Fatal error when installing Python module using PIP.\")
     endif()"
  )
  # Copy the setup file and README.rst into the binary dir so you can install
  # using pip or make a distribution with the compiled binary
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Python/pkginit.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_kinetx
  )
  file(RENAME
    ${CMAKE_CURRENT_BINARY_DIR}/scine_kinetx/pkginit.py
    ${CMAKE_CURRENT_BINARY_DIR}/scine_kinetx/__init__.py
  )
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Python/README.rst DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  # Typing stubs
  include(FindPythonModule)
  find_python_module(pybind11_stubgen)
  if(PY_PYBIND11_STUBGEN)
    add_custom_command(TARGET scine_kinetx POST_BUILD
      COMMAND ${PYTHON_EXECUTABLE} -m pybind11_stubgen -o . --root-module-suffix \"\" --no-setup-py --bare-numpy-ndarray scine_kinetx
      COMMENT "Generating python package typechecking stubs with pybind11-stubgen"
    )
  else()
    message(STATUS "Not generating typechecking stubs for python package as pybind11-stubgen was not found")
  endif()

  # Add setuptools file
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
    ${CMAKE_CURRENT_BINARY_DIR}/setup.py
    @ONLY
  )
  file(
    COPY ${PROJECT_SOURCE_DIR}/README.rst
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )

  # Python Binding Tests
  if(SCINE_BUILD_TESTS)
    add_test(
      NAME Kinetx_Python_Bindings
      COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/Python/Tests --junitxml=${CMAKE_CURRENT_BINARY_DIR}/pytest_report.xml
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(Kinetx_Python_Bindings PROPERTIES
      ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()

  include(SphinxDocumentation)
  scine_sphinx_documentation(
    TARGET scine_kinetx
    CONFIGURATION ${CMAKE_CURRENT_SOURCE_DIR}/Python/Sphinx/conf.py.in
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Python/Sphinx
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doc-py
    LINK UtilsPythonDoc
    DOCTEST
  )

  file(COPY ${PROJECT_SOURCE_DIR}/LICENSE.txt
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_kinetx/
  )
  file(COPY ${EXTERNAL_LICENSE_FILES}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_kinetx/
  )
endif()
