// Python-Modul cliffordbem._cbem (v0.42): Einstiegspunkt, Version, Threads.
#include "common.hpp"

#ifndef CBEM_VERSION
#define CBEM_VERSION "unbekannt"
#endif

PYBIND11_MODULE(_cbem, m) {
    namespace py = pybind11;
    using namespace cbem::py_bind;
    m.doc() = "CliffordBem: koordinatenfreie Randelementmethode der Nano-Optik in Cl3(C) (C++-Kern)";
    m.attr("__version__") = CBEM_VERSION;
#ifdef CBEM_USE_OPENMP
    m.attr("HAS_OPENMP") = true;
#else
    m.attr("HAS_OPENMP") = false;
#endif
    m.def("omp_threads", &cbem::omp_threads, "Zahl der OpenMP-Threads (1 ohne OpenMP)");
    m.def("set_num_threads", [](int n) {
              if (n < 1) throw py::value_error("mindestens ein Thread");
#ifdef CBEM_USE_OPENMP
              omp_set_num_threads(n);
#endif
          }, "n"_a, "Zahl der OpenMP-Threads setzen (ohne OpenMP ohne Wirkung)");
    m.attr("pi") = cbem::pi;
#ifdef CBEM_SOURCE_DATA_DIR
    m.attr("_SOURCE_DATA_DIR") = CBEM_SOURCE_DATA_DIR;       // Materialdaten im Quellbaum (Bau mit CMake)
#else
    m.attr("_SOURCE_DATA_DIR") = "";
#endif

    // Reihenfolge: Grundtypen (Medium, Parameter) vor den Klassen, die sie als Voreinstellung verwenden
    init_core(m);
    init_geometry(m);
    init_operators(m);
    init_problems(m);
    init_sources(m);
    init_linear(m);
    init_curved(m);
}
