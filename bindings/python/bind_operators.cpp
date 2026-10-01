// Python-Anbindung: Eintragsauswertung, H-Matrix, Randoperatoren, Transmissionsoperator, GMRES (v0.42).
//
// Diese Klassen halten intern Referenzen auf ihre Eingaben (Netz, Eintraege, H-Matrix, Operatoren). Die Anbindung haelt die
// Eingaben deshalb ueber keep_alive am Leben, solange das abhaengige Objekt existiert.
#include "common.hpp"

#include "cbem/assembly/kernel_entries.hpp"
#include "cbem/hmatrix/hmatrix.hpp"
#include "cbem/operators/cauchy_operator.hpp"
#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/solvers/gmres.hpp"

namespace cbem::py_bind {

namespace {

py::array_t<cplx> comp_to_numpy(const KernelComp& K) {
    py::array_t<cplx> a(4); std::copy(K.begin(), K.end(), a.mutable_data()); return a;
}

// Python-Funktion x -> y als LinOp; holt sich bei jedem Aufruf den GIL (GMRES selbst laeuft ohne GIL)
LinOp python_linop(const py::function& f, std::size_t n, const char* what) {
    // shared_ptr: Kopien des LinOp ohne GIL beruehren den Python-Referenzzaehler nicht
    auto fp = std::make_shared<py::function>(f);
    return [fp, n, what](const std::vector<cplx>& x, std::vector<cplx>& y) {
        py::gil_scoped_acquire gil;
        py::object r = (*fp)(to_numpy(x));
        CArr a = r.cast<CArr>();
        if (static_cast<std::size_t>(a.size()) != n)
            throw py::value_error(std::string(what) + ": Ergebnis hat Laenge " + std::to_string(a.size()) + " statt " + std::to_string(n));
        y.assign(a.data(), a.data() + a.size());
    };
}

}  // namespace

void init_operators(py::module_& m) {
    py::class_<KernelEntries>(m, "KernelEntries", R"doc(
Eintragsauswertung des Cauchy-Operators E_k: Kernkomponenten K(i, j) = int_tau_i int_tau_j Phi_k(x - y) dS_y dS_x
(Skalar + Vektor, 4 komplexe Werte). Haelt das Netz am Leben.
)doc")
        .def(py::init<const TriangleMesh&, cplx, EntryParams>(), "mesh"_a, "k"_a, "params"_a = EntryParams{}, py::keep_alive<1, 2>(),
             py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("wavenumber", &KernelEntries::wavenumber)
        .def("exact", [](const KernelEntries& E, std::size_t i, std::size_t j) { return comp_to_numpy(E.exact(i, j)); }, "i"_a, "j"_a,
             "genaue Eintraege (Nahfeld, Sauter-Schwab oder Gauss je nach Paar)")
        .def("far", [](const KernelEntries& E, std::size_t i, std::size_t j) { return comp_to_numpy(E.far(i, j)); }, "i"_a, "j"_a,
             "Gauss-Fernfeldregel (nur fuer getrennte Paare genau)")
        .def("near", [](const KernelEntries& E, std::size_t i, std::size_t j) { return comp_to_numpy(E.near(i, j)); }, "i"_a, "j"_a)
        .def("is_near", &KernelEntries::is_near, "i"_a, "j"_a)
        .def("adjacency", &KernelEntries::adjacency, "i"_a, "j"_a)
        .def("aspect", &KernelEntries::aspect, "t"_a, "Seitenverhaeltnis h^2 / (2 A)")
        .def_property_readonly("near_pairs", &KernelEntries::near_pairs)
        .def_property_readonly("near_seconds", &KernelEntries::near_seconds);

    py::class_<KernelHMatrix>(m, "KernelHMatrix", "H-Matrix der vier Kernkomponenten (ACA)")
        .def(py::init<const KernelEntries&, HMatrixParams>(), "entries"_a, "params"_a = HMatrixParams{}, py::keep_alive<1, 2>(),
             py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("stats", &KernelHMatrix::stats)
        .def("__len__", &KernelHMatrix::size)
        .def("apply", [](const KernelHMatrix& H, const CArr& Z) {
                 if (static_cast<std::size_t>(Z.size()) != 32 * H.size()) throw py::value_error("Z: Laenge 32 N erwartet (N x 4 x 8)");
                 auto z = to_cvec(Z);
                 std::vector<cplx> y(8 * H.size(), cplx(0));
                 nogil([&] { H.apply(z, y); });
                 return to_numpy(y);
             }, "Z"_a, "Y_i = sum_c sum_j K_c(i, j) Z_(j, c), Z: N x 4 x 8, Y: N x 8");

    py::class_<BoundaryOperator>(m, "BoundaryOperator", "Randoperator auf der stueckweise konstanten Basis (8 Komponenten je Dreieck)")
        .def("apply", [](const BoundaryOperator& E, const CArr& x) {
                 auto xv = to_trace(x, E.size() / 8, "x");
                 std::vector<cplx> y;
                 nogil([&] { E.apply(xv, y); });
                 return to_numpy(y);
             }, "x"_a)
        .def("__matmul__", [](const BoundaryOperator& E, const CArr& x) {
            auto xv = to_trace(x, E.size() / 8, "x");
            std::vector<cplx> y;
            nogil([&] { E.apply(xv, y); });
            return to_numpy(y);
        })
        .def("__len__", &BoundaryOperator::size)
        .def_property_readonly("shape", [](const BoundaryOperator& E) { return py::make_tuple(E.size(), E.size()); });
    py::class_<CauchyOperator, BoundaryOperator>(m, "CauchyOperator", "Cauchy-Randoperator E_k (Involution: E E = 1)")
        .def(py::init<const TriangleMesh&, const KernelHMatrix&>(), "mesh"_a, "H"_a, py::keep_alive<1, 2>(), py::keep_alive<1, 3>());
    py::class_<ChiralCauchyOperator, BoundaryOperator>(m, "ChiralCauchyOperator", "chiraler Innenoperator P+ E_(k+) + P- E_(k-)")
        .def(py::init<const CauchyOperator&, const CauchyOperator&>(), "E_plus"_a, "E_minus"_a, py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

    py::class_<TransmissionOperator>(m, "TransmissionOperator", "T_1 = 1/2 (1 + E_2) + 1/2 (1 - E_1) J")
        .def(py::init<const TriangleMesh&, const BoundaryOperator&, const BoundaryOperator&, const Medium&, const Medium&>(), "mesh"_a,
             "E_inner"_a, "E_outer"_a, "inner"_a, "outer"_a, py::keep_alive<1, 3>(), py::keep_alive<1, 4>())
        .def(py::init<const TriangleMesh&, const BoundaryOperator&, const BoundaryOperator&, const std::vector<Medium>&, const Medium&>(),
             "mesh"_a, "E_inner"_a, "E_outer"_a, "inner_per_triangle"_a, "outer"_a, py::keep_alive<1, 3>(), py::keep_alive<1, 4>())
        .def("apply", [](const TransmissionOperator& T, const CArr& x) {
                 auto xv = to_trace(x, T.size() / 8, "x");
                 std::vector<cplx> y;
                 nogil([&] { T.apply(xv, y); });
                 return to_numpy(y);
             }, "x"_a)
        .def("precondition", [](const TransmissionOperator& T, const CArr& x) {
                 auto xv = to_trace(x, T.size() / 8, "x");
                 std::vector<cplx> y;
                 nogil([&] { T.precondition(xv, y); });
                 return to_numpy(y);
             }, "x"_a, "punktweise Vorkonditionierung 2 (1 + J)^-1")
        .def("__len__", &TransmissionOperator::size)
        .def_property_readonly("J", [](const TransmissionOperator& T) {
            const auto& J = T.J();
            py::array_t<cplx> a({static_cast<py::ssize_t>(J.size()), py::ssize_t(8), py::ssize_t(8)});
            for (std::size_t t = 0; t < J.size(); ++t) std::copy(J[t].begin(), J[t].end(), a.mutable_data() + 64 * t);
            return a;
        }, "Transmissionsabbildungen je Dreieck (N, 8, 8)");

    py::class_<GmresResult>(m, "GmresResult")
        .def_readonly("iterations", &GmresResult::iterations)
        .def_readonly("rel_residual", &GmresResult::rel_residual)
        .def_readonly("converged", &GmresResult::converged)
        .def("__repr__", [](const GmresResult& r) {
            return "GmresResult(iterations=" + std::to_string(r.iterations) + ", rel_residual=" +
                   py::repr(py::float_(r.rel_residual)).cast<std::string>() + ", converged=" + (r.converged ? "True" : "False") + ")";
        });
    m.def("gmres", [](py::object A, const CArr& b, py::object M, real tol, int restart, int max_iter) {
              const auto bv = to_cvec(b);
              const std::size_t n = bv.size();
              // C++-Operatoren direkt (ohne Umweg ueber Python), sonst beliebige Funktion x -> y
              LinOp Aop, Mop;
              if (py::isinstance<TransmissionOperator>(A)) {
                  const auto* T = A.cast<const TransmissionOperator*>();
                  Aop = [T](const std::vector<cplx>& x, std::vector<cplx>& y) { T->apply(x, y); };
              } else if (py::isinstance<BoundaryOperator>(A)) {
                  const auto* E = A.cast<const BoundaryOperator*>();
                  Aop = [E](const std::vector<cplx>& x, std::vector<cplx>& y) { E->apply(x, y); };
              } else Aop = python_linop(A.cast<py::function>(), n, "A");
              if (!M.is_none()) Mop = python_linop(M.cast<py::function>(), n, "M");
              std::vector<cplx> x;
              const LinOp* Mp = M.is_none() ? nullptr : &Mop;
              GmresResult r = nogil([&] { return gmres(Aop, bv, x, Mp, tol, restart, max_iter); });
              return py::make_tuple(to_numpy(x), r);
          }, "A"_a, "b"_a, "M"_a = py::none(), "tol"_a = 1e-6, "restart"_a = 200, "max_iter"_a = 2000, R"doc(
Neugestartetes GMRES mit Rechtsvorkonditionierung: loest A x = b. A ist ein TransmissionOperator, ein
BoundaryOperator oder eine Funktion x -> A x (NumPy-Arrays); M (optional) eine Funktion x -> M x.
Gibt (x, GmresResult) zurueck.
)doc");
}

}  // namespace cbem::py_bind
