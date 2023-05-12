cmake_minimum_required(VERSION 3.9)

project(KinetxTestPackage)

find_package(ScineKinetx REQUIRED)

add_executable(TestPackageTest ${CMAKE_CURRENT_SOURCE_DIR}/test.cpp)
target_link_libraries(TestPackageTest PRIVATE Scine::Kinetx)
enable_testing()
add_test(NAME LinkingTest COMMAND TestPackageTest)
