set(KINETX_APP_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/App/main.cpp
)

set(KINETX_LIB_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/Network.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/Network.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/NetworkBuilder.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/NetworkBuilder.h
	${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/RandomNetworkFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/RandomNetworkFactory.h
	${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/ReferenceNetworks.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/ReferenceNetworks.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/RungeKutta/CashKarp5.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/RungeKutta/CashKarp5.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/RungeKutta/ExplicitEuler.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/RungeKutta/ExplicitEuler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/RungeKutta/ImplicitEuler.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/RungeKutta/ImplicitEuler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/RungeKutta/RungeKutta.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Kinetx/RungeKutta/RungeKutta.cpp
)

set(KINETX_TEST_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/Tests/NetworkBuilderTest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Tests/NetworkTest.cpp
)

set(KINETX_PYTHON_CPPS
  ${CMAKE_CURRENT_SOURCE_DIR}/Python/PythonModule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Python/NetworkPython.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Python/NetworkBuilderPython.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Python/ReferenceNetworksPython.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Python/RandomNetworkFactoryPython.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Python/RungeKuttaPython.cpp
)
