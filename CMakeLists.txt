cmake_minimum_required(VERSION 3.16)
project(hlsimii CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(hlsim_core STATIC src/line.cpp src/mesh.cpp src/config.cpp src/physics.cpp src/linalg.cpp src/fall_simulator.cpp)
target_include_directories(hlsim_core PUBLIC src)
target_link_libraries(hlsim_core PUBLIC Eigen3::Eigen)

# YAML front end for SimConfig; the WASM build links only hlsim_core
add_library(hlsim_yaml STATIC src/config_yaml.cpp)
target_link_libraries(hlsim_yaml PUBLIC hlsim_core)

add_executable(hlsimii src/main.cpp)
target_link_libraries(hlsimii PRIVATE hlsim_yaml)

enable_testing()
add_executable(test_line_mesh tests/test_line_mesh.cpp)
target_link_libraries(test_line_mesh PRIVATE hlsim_core)
add_test(NAME line_mesh COMMAND test_line_mesh)

add_executable(test_config tests/test_config.cpp)
target_link_libraries(test_config PRIVATE hlsim_yaml)
target_compile_definitions(test_config PRIVATE HLSIM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
add_test(NAME config COMMAND test_config)

add_executable(test_physics tests/test_physics.cpp)
target_link_libraries(test_physics PRIVATE hlsim_core)
add_test(NAME physics COMMAND test_physics)

add_executable(test_linalg tests/test_linalg.cpp)
target_link_libraries(test_linalg PRIVATE hlsim_core)
add_test(NAME linalg COMMAND test_linalg)

add_executable(test_free_fall tests/test_free_fall.cpp)
target_link_libraries(test_free_fall PRIVATE hlsim_core)
add_test(NAME free_fall COMMAND test_free_fall)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  foreach(target hlsim_core hlsim_yaml hlsimii test_line_mesh test_config test_physics test_linalg test_free_fall)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endforeach()
endif()
