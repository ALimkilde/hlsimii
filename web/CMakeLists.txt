cmake_minimum_required(VERSION 3.16)
project(hlsimii_web CXX)

# Builds the hlsimii simulator to WebAssembly and assembles the static site in dist/.
# This is a separate project from the native build in ../CMakeLists.txt (it needs emcmake).
# It lists the simulator sources directly instead of add_subdirectory(..), which would also
# build main.cpp and the tests and look for Eigen with find_package.

if(NOT EMSCRIPTEN)
  message(FATAL_ERROR "Configure with emcmake: emcmake cmake -B build")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(HLSIM_DIR "${CMAKE_SOURCE_DIR}/.." CACHE PATH "hlsimii source root")
# Eigen is header-only; find_package does not look outside the emscripten sysroot
set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3" CACHE PATH "Eigen headers")
# worker: for the web page; node: for test/smoke.mjs
set(HLSIM_ENVIRONMENT "worker" CACHE STRING "emscripten -sENVIRONMENT")
set(HLSIM_OUT_DIR "${CMAKE_SOURCE_DIR}/dist" CACHE PATH "where hlsim.mjs/.wasm go: the published site")

add_executable(hlsim
  ${HLSIM_DIR}/src/line.cpp
  ${HLSIM_DIR}/src/mesh.cpp
  ${HLSIM_DIR}/src/config.cpp
  ${HLSIM_DIR}/src/config_yaml.cpp
  ${HLSIM_DIR}/src/physics.cpp
  ${HLSIM_DIR}/src/linalg.cpp
  ${HLSIM_DIR}/src/fall_simulator.cpp
  ${HLSIM_DIR}/src/experiments.cpp
  bindings.cpp)
target_include_directories(hlsim PRIVATE ${HLSIM_DIR}/src ${EIGEN3_INCLUDE_DIR})
# validate() and the YAML parser report errors by throwing; wasm exceptions keep them catchable
target_compile_options(hlsim PRIVATE -O3 -fwasm-exceptions)
target_link_options(hlsim PRIVATE
  -O3 -lembind -fwasm-exceptions
  -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createHlsim
  -sENVIRONMENT=${HLSIM_ENVIRONMENT}
  -sALLOW_MEMORY_GROWTH=1)
set_target_properties(hlsim PROPERTIES
  SUFFIX ".mjs"
  RUNTIME_OUTPUT_DIRECTORY ${HLSIM_OUT_DIR})

# dist/ = the page in site/ + the simulator's webbing catalogue and examples + the documentation PDF
# (the page's Documentation button links to it). A target of its
# own so edits there reach dist/ on every build, not only when hlsim relinks. Skipped for the
# node build (test/smoke.mjs reads the YAML from the source tree). The examples must match the
# <select id="example"> options in site/index.html.
set(HLSIM_EXAMPLES 110m.yaml short_and_static.yaml 250m.yaml)
if(HLSIM_ENVIRONMENT STREQUAL "worker")
  list(TRANSFORM HLSIM_EXAMPLES PREPEND ${HLSIM_DIR}/examples/ OUTPUT_VARIABLE example_paths)
  add_custom_target(site ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/site ${HLSIM_OUT_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy ${HLSIM_DIR}/data/webbings.yaml ${example_paths} ${HLSIM_OUT_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy ${HLSIM_DIR}/docs/hlsim_numerical_model.pdf ${HLSIM_OUT_DIR}
    VERBATIM)
endif()
