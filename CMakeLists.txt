cmake_minimum_required(VERSION 3.16)
project(CliffordBem VERSION 0.32.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(CBEM_OPENMP "OpenMP-Parallelisierung der Assemblierung" ON)
option(CBEM_NATIVE "Mit -march=native uebersetzen" OFF)
option(CBEM_BUILD_TESTS "Tests bauen" ON)
option(CBEM_BUILD_APPS "Anwendungen/Benchmarks bauen" ON)

add_library(cbem
  src/multivector.cpp
  src/mesh.cpp
  src/gmsh_io.cpp
  src/quadrature.cpp
  src/sauter_schwab.cpp
  src/triangle_integrals.cpp
  src/kernel_entries.cpp
  src/dense.cpp
  src/cluster_tree.cpp
  src/aca.cpp
  src/hmatrix.cpp
  src/cauchy_operator.cpp
  src/chiral_cauchy_operator.cpp
  src/transmission_operator.cpp
  src/gmres.cpp
  src/fields.cpp
  src/dense_blocks.cpp
  src/block_preconditioner.cpp
  src/scattering_problem.cpp
  src/materials.cpp
  src/orientation.cpp
  src/hodlr.cpp
  src/layered_problem.cpp
  src/thin_layer_problem.cpp
  src/twoport_layer_problem.cpp
  src/chiral_incidence.cpp
  src/near_field.cpp
  src/dipole.cpp
  src/recycling_gmres.cpp)
target_include_directories(cbem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(cbem PRIVATE -Wall -Wextra)
if(CBEM_NATIVE)
  target_compile_options(cbem PUBLIC -march=native)
endif()
if(CBEM_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(cbem PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(cbem PUBLIC CBEM_USE_OPENMP)
  endif()
endif()

if(CBEM_BUILD_APPS)
  add_executable(bench_compression apps/bench_compression.cpp)
  target_link_libraries(bench_compression PRIVATE cbem)
  add_executable(scatter_sphere apps/scatter_sphere.cpp)
  target_link_libraries(scatter_sphere PRIVATE cbem)
  add_executable(scatter_cube apps/scatter_cube.cpp)
  target_link_libraries(scatter_cube PRIVATE cbem)
  add_executable(scatter_chiral apps/scatter_chiral.cpp)
  target_link_libraries(scatter_chiral PRIVATE cbem)
  add_executable(scatter_multi apps/scatter_multi.cpp)
  target_link_libraries(scatter_multi PRIVATE cbem)
  add_executable(scatter_mesh apps/scatter_mesh.cpp)
  target_link_libraries(scatter_mesh PRIVATE cbem)
  add_executable(dipole apps/dipole.cpp)
  target_link_libraries(dipole PRIVATE cbem)
  add_executable(nearfield apps/nearfield.cpp)
  target_link_libraries(nearfield PRIVATE cbem)
  add_executable(spectrum apps/spectrum.cpp)
  target_link_libraries(spectrum PRIVATE cbem)
  add_executable(scatter_cube_dimer apps/scatter_cube_dimer.cpp)
  target_link_libraries(scatter_cube_dimer PRIVATE cbem)
  add_executable(scatter_coated apps/scatter_coated.cpp)
  target_link_libraries(scatter_coated PRIVATE cbem)
endif()

if(CBEM_BUILD_TESTS)
  enable_testing()
  foreach(t test_multivector test_triangle_integrals test_dense test_aca test_hmatrix test_plemelj test_gmres test_scattering test_block_precond test_sauter_schwab test_nearfield test_chiral test_multibody test_gmsh test_materials test_hodlr test_layered test_thin_layer test_twoport test_chiral_host test_near_field test_dipole test_recycling)
    add_executable(${t} tests/${t}.cpp)
    target_link_libraries(${t} PRIVATE cbem)
    add_test(NAME ${t} COMMAND ${t})
  endforeach()
  target_compile_definitions(test_materials PRIVATE CBEM_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/materials")
endif()
