cmake_minimum_required(VERSION 3.16)
project(CliffordBem VERSION 0.54.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(CBEM_OPENMP "OpenMP-Parallelisierung der Assemblierung" ON)
option(CBEM_NATIVE "Mit -march=native uebersetzen" OFF)
option(CBEM_BUILD_TESTS "Tests bauen" ON)
option(CBEM_BUILD_APPS "Anwendungen/Benchmarks bauen" ON)
option(CBEM_BUILD_PYTHON "Python-Modul cliffordbem (pybind11) bauen" OFF)

add_library(cbem
  src/multivector.cpp
  src/mesh.cpp
  src/gmsh_io.cpp
  src/quadrature.cpp
  src/sauter_schwab.cpp
  src/triangle_integrals.cpp
  src/kernel_entries.cpp
  src/linear_entries.cpp
  src/linear_cauchy_operator.cpp
  src/linear_problem.cpp
  src/quadratic_mesh.cpp
  src/curved_entries.cpp
  src/curved_hmatrix.cpp
  src/curved_operators.cpp
  src/curved_problem.cpp
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
  src/recycling_gmres.cpp
  src/optical_force.cpp
  src/incident_field.cpp)
target_include_directories(cbem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(cbem PRIVATE -Wall -Wextra)
if(CBEM_BUILD_PYTHON)
  set_target_properties(cbem PROPERTIES POSITION_INDEPENDENT_CODE ON)   # wird in das Python-Modul gelinkt
endif()
if(CBEM_NATIVE)
  target_compile_options(cbem PUBLIC -march=native)
endif()
# MinGW (Windows): Laufzeitbibliotheken (libstdc++, libgcc, libgomp, winpthread) statisch linken; sonst laden die
# Programme je nach PATH fremde, inkompatible DLLs (z. B. libstdc++-6.dll aus Git for Windows) und stuerzen ab.
# Statische Archive bevorzugen, damit find_package(OpenMP) libgomp.a statt libgomp.dll.a eintraegt.
if(MINGW)
  add_link_options(-static)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
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
  add_executable(tweezers apps/tweezers.cpp)
  target_link_libraries(tweezers PRIVATE cbem)
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
  foreach(t test_multivector test_triangle_integrals test_dense test_aca test_hmatrix test_plemelj test_gmres test_scattering test_block_precond test_sauter_schwab test_nearfield test_chiral test_multibody test_gmsh test_materials test_hodlr test_layered test_thin_layer test_twoport test_chiral_host test_near_field test_dipole test_recycling test_optical_force test_beams test_linear test_curved)
    add_executable(${t} tests/${t}.cpp)
    target_link_libraries(${t} PRIVATE cbem)
    add_test(NAME ${t} COMMAND ${t})
  endforeach()
  target_compile_definitions(test_materials PRIVATE CBEM_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/materials")
  target_compile_definitions(test_gmsh PRIVATE CBEM_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/meshes")
  # GLMT-Referenz (v0.39): Selbstpruefung des Python-Werkzeugs, falls Python 3 mit numpy und scipy vorhanden ist
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_FOUND)
    add_test(NAME test_glmt COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/glmt.py --selftest)
  endif()
endif()

# Python-Anbindung (v0.42): Modul cliffordbem._cbem. Benoetigt pybind11 (pip install pybind11) und die Python-Entwicklungsdateien.
# Nach dem Bau liegt das Paket in build/python (PYTHONPATH=build/python); pip install . baut es ueber pyproject.toml.
if(CBEM_BUILD_PYTHON)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  if(NOT DEFINED pybind11_DIR)
    execute_process(COMMAND ${Python_EXECUTABLE} -m pybind11 --cmakedir OUTPUT_VARIABLE _cbem_pybind11_dir
                    OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE _cbem_pybind11_res ERROR_QUIET)
    if(_cbem_pybind11_res EQUAL 0)
      list(APPEND CMAKE_PREFIX_PATH "${_cbem_pybind11_dir}")
    endif()
  endif()
  set(PYBIND11_FINDPYTHON ON)
  find_package(pybind11 2.10 CONFIG REQUIRED)
  pybind11_add_module(_cbem MODULE
    bindings/python/module.cpp
    bindings/python/bind_core.cpp
    bindings/python/bind_geometry.cpp
    bindings/python/bind_operators.cpp
    bindings/python/bind_problems.cpp
    bindings/python/bind_sources.cpp
    bindings/python/bind_linear.cpp
    bindings/python/bind_curved.cpp)
  target_link_libraries(_cbem PRIVATE cbem)
  target_compile_options(_cbem PRIVATE -Wall -Wextra)
  target_compile_definitions(_cbem PRIVATE CBEM_VERSION="${PROJECT_VERSION}"
                                           CBEM_SOURCE_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/materials")
  set(CBEM_PY_DIR ${CMAKE_BINARY_DIR}/python/cliffordbem)
  set_target_properties(_cbem PROPERTIES LIBRARY_OUTPUT_DIRECTORY "$<1:${CBEM_PY_DIR}>")   # ohne Unterverzeichnis je Konfiguration
  file(GLOB _cbem_py_files CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/bindings/python/cliffordbem/*.py)
  foreach(_f ${_cbem_py_files})
    get_filename_component(_n ${_f} NAME)
    configure_file(${_f} ${CBEM_PY_DIR}/${_n} COPYONLY)
  endforeach()
  if(SKBUILD)   # pip install . (scikit-build-core): Modul und Materialdaten ins Paket
    install(TARGETS _cbem LIBRARY DESTINATION cliffordbem)
    install(DIRECTORY data/materials DESTINATION cliffordbem/data)
  endif()
  if(CBEM_BUILD_TESTS)
    add_test(NAME test_python COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/bindings/python/tests/test_python.py)
    # OMP_NUM_THREADS=1: die Tests vergleichen Wege bitgenau (np.array_equal); mit mehreren Threads summiert die H-Matrix in
    # nicht festgelegter Reihenfolge
    set_tests_properties(test_python PROPERTIES ENVIRONMENT
      "PYTHONPATH=${CMAKE_BINARY_DIR}/python;CBEM_TOOLS_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tools;OMP_NUM_THREADS=1")
  endif()
endif()
