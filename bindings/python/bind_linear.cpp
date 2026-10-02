// Python-Anbindung: unstetig lineare Dichten auf ebenen Dreiecken (v0.45, Stufe 2a der gekruemmten Elemente).
// Spuren haben 24 Koeffizienten je Dreieck: h[8 (3 t + a) + q], a = 0..2 die je Element orthonormierte Basis psi_a,
// q das Blade. reshape(-1, 3, 8) gibt (Dreieck, Basisfunktion, Blade).
#include "common.hpp"

#include "cbem/problems/linear_problem.hpp"
#include "cbem/sources/chiral_incidence.hpp"

namespace cbem::py_bind {

namespace {

py::array_t<cplx> block_to_numpy(const LinearBlock& B) {
    py::array_t<cplx> a({py::ssize_t(3), py::ssize_t(3), py::ssize_t(4)});
    cplx* d = a.mutable_data();
    for (int p = 0; p < 9; ++p) for (int c = 0; c < 4; ++c) d[4 * p + c] = B[p][c];
    return a;
}

std::vector<cplx> to_linear_trace(const CArr& a, std::size_t N, const char* what) {
    if (static_cast<std::size_t>(a.size()) != 24 * N)
        throw py::value_error(std::string(what) + ": erwartet 24 N = " + std::to_string(24 * N) + " Koeffizienten (lineare Dichten), erhalten " +
                              std::to_string(a.size()));
    return to_cvec(a);
}

void check_media(std::size_t bodies, std::size_t media) {
    if (bodies == 0) throw py::value_error("mindestens ein Koerper erwartet");
    if (bodies != media) throw py::value_error("je Koerper ein Medium erwartet");
}

}  // namespace

void init_linear(py::module_& m) {
    py::class_<LinearKernelEntries>(m, "LinearKernelEntries", R"doc(
Eintragsauswertung des Cauchy-Operators fuer unstetig lineare Dichten: je Elementpaar 3 x 3 Formfunktionen x 4
Kernkomponenten (Skalar + Vektor). block(i, j) in der je Element orthonormierten Basis psi, lambda_block(i, j) in den
baryzentrischen Koordinaten. Summe ueber die Formfunktionen = KernelEntries. Haelt das Netz am Leben.
)doc")
        .def(py::init<const TriangleMesh&, cplx, EntryParams>(), "mesh"_a, "k"_a, "params"_a = EntryParams{}, py::keep_alive<1, 2>(),
             py::call_guard<py::gil_scoped_release>())
        .def("block", [](const LinearKernelEntries& E, std::size_t i, std::size_t j) { return block_to_numpy(E.block(i, j)); }, "i"_a, "j"_a,
             "psi-Basis, Form (3, 3, 4): [a, b, c]")
        .def("lambda_block", [](const LinearKernelEntries& E, std::size_t i, std::size_t j) { return block_to_numpy(E.lambda_exact(i, j)); },
             "i"_a, "j"_a, "baryzentrische Basis lambda, Form (3, 3, 4)")
        .def("S", [](const LinearKernelEntries& E, std::size_t t) {
                 py::array_t<double> a({py::ssize_t(3), py::ssize_t(3)}); std::copy(E.S(t).begin(), E.S(t).end(), a.mutable_data()); return a;
             }, "t"_a, "psi_a = sum_k S[a, k] lambda_k")
        .def_property_readonly("near_pairs", &LinearKernelEntries::near_pairs)
        .def("__len__", &LinearKernelEntries::size);

    py::class_<LinearCauchyOperator, BoundaryOperator>(m, "LinearCauchyOperator", "Cauchy-Operator fuer lineare Dichten (24 N)")
        .def(py::init<const TriangleMesh&, const KernelHMatrix&>(), "mesh"_a, "H"_a, py::keep_alive<1, 2>(), py::keep_alive<1, 3>());
    // H-Matrix aus LinearKernelEntries (KernelHMatrix ist bereits gebunden; zusaetzlicher Konstruktor)
    m.def("linear_hmatrix", [](const LinearKernelEntries& E, const HMatrixParams& hp) {
              return nogil([&] { return std::make_unique<KernelHMatrix>(E, hp); });
          }, "entries"_a, "params"_a = HMatrixParams{}, py::keep_alive<0, 1>(), "KernelHMatrix mit 3 N Indizes (nur AcaMode.Joint)");

    m.def("project_plane_wave_linear", [](const TriangleMesh& me, cplx k, cplx eps, const Vec3& d, const CVec3& p, int sub) {
              return to_numpy(nogil([&] { return project_plane_wave_linear(me, k, eps, d, p, sub); }));
          }, "mesh"_a, "k"_a, "eps"_a, "d"_a, "p"_a, "sub"_a = 2, "L2-Projektion der ebenen Welle auf die linearen Dichten (24 N)");
    m.def("plane_wave_trace_linear", [](const TriangleMesh& me, const Medium& md, real omega, const Vec3& d, const CVec3& p, int sub) {
              const cplx k = plane_wave_incidence(md, omega, d, p).k;     // chiral: Helizitaetswelle mit k_sigma (v0.59)
              return to_numpy(nogil([&] { return project_plane_wave_linear(me, k, md.eps, d, p, sub); }));
          }, "mesh"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "sub"_a = 2, "rechte Seite wie LinearScatteringProblem.solve_plane_wave");
    m.def("far_field_linear", [](const TriangleMesh& me, const CArr& hs, cplx k, const Vec3& xhat, int sub) {
              auto h = to_linear_trace(hs, me.size(), "h_scat");
              return nogil([&] { return far_field_linear(me, h, k, xhat, sub); });
          }, "mesh"_a, "h_scat"_a, "k"_a, "xhat"_a, "sub"_a = 2);
    m.def("extinction_cross_section_linear", [](const TriangleMesh& me, const CArr& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
              auto h = to_linear_trace(hs, me.size(), "h_scat");
              return nogil([&] { return extinction_cross_section_linear(me, h, k, eps, d, p); });
          }, "mesh"_a, "h_scat"_a, "k"_a, "eps"_a, "d"_a, "p"_a);
    m.def("forward_amplitude_linear", [](const TriangleMesh& me, const CArr& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
              auto h = to_linear_trace(hs, me.size(), "h_scat");
              return nogil([&] { return forward_amplitude_linear(me, h, k, eps, d, p); });
          }, "mesh"_a, "h_scat"_a, "k"_a, "eps"_a, "d"_a, "p"_a);
    m.def("linear_to_constant", [](const TriangleMesh& me, const CArr& h) { return to_numpy(linear_to_constant(me, to_linear_trace(h, me.size(), "h"))); },
          "mesh"_a, "h"_a, "Mittelwert je Dreieck als konstante Spur (8 N, Normierung wie ScatteringProblem)");
    m.def("linear_trace_value", [](const TriangleMesh& me, const CArr& h, std::size_t t, const Vec3& x) {
              if (t >= me.size()) throw py::index_error("Dreiecksindex");
              return linear_trace_value(me, to_linear_trace(h, me.size(), "h"), t, x);
          }, "mesh"_a, "h"_a, "t"_a, "x"_a, "Dichte (Multivektor) im Punkt x des Dreiecks t");

    py::class_<LinearScatteringProblem>(m, "LinearScatteringProblem", R"doc(
Streuproblem mit unstetig linearen Dichten auf ebenen Dreiecken (24 Unbekannte je Dreieck), sonst wie ScatteringProblem:
ein oder mehrere Koerper (auch chiral), achirales Aussenmedium, ebene Wellen und beliebige rechte Seiten (Laenge 24 N).
Auf ebenen Elementen allein bringt das an glatten Koerpern keinen Gewinn (der Geometriefehler dominiert,
docs/results_curved.md); es ist die Grundlage der gekruemmten Elemente.
)doc")
        .def(py::init([](const std::vector<TriangleMesh>& bodies, const std::vector<Medium>& media, real omega, const Medium& outer,
                         const HMatrixParams& hp, const EntryParams& ep) {
                 check_media(bodies.size(), media.size());
                 return nogil([&] { return std::make_unique<LinearScatteringProblem>(bodies, media, omega, outer, hp, ep); });
             }),
             "bodies"_a, "media"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = HMatrixParams{}, "entries"_a = EntryParams{})
        .def(py::init([](const TriangleMesh& body, const Medium& medium, real omega, const Medium& outer, const HMatrixParams& hp,
                         const EntryParams& ep) {
                 return nogil([&] {
                     return std::make_unique<LinearScatteringProblem>(std::vector<TriangleMesh>{body}, std::vector<Medium>{medium}, omega, outer, hp, ep);
                 });
             }),
             "body"_a, "medium"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = HMatrixParams{}, "entries"_a = EntryParams{})
        .def("solve_plane_wave", &LinearScatteringProblem::solve_plane_wave, "d"_a, "p"_a, "options"_a = SolveOptions{},
             py::call_guard<py::gil_scoped_release>())
        .def("solve_rhs", [](const LinearScatteringProblem& P, const CArr& b, const SolveOptions& o) {
                 auto bv = to_linear_trace(b, P.mesh().size(), "b");
                 return nogil([&] { return P.solve_rhs(bv, o); });
             }, "b"_a, "options"_a = SolveOptions{})
        .def("apply", [](const LinearScatteringProblem& P, const CArr& x) {
                 auto xv = to_linear_trace(x, P.mesh().size(), "x");
                 std::vector<cplx> y; nogil([&] { P.T().apply(xv, y); }); return to_numpy(y);
             }, "x"_a)
        .def_property_readonly("mesh", &LinearScatteringProblem::mesh, py::return_value_policy::reference_internal)
        .def_property_readonly("body_begin", [](const LinearScatteringProblem& P) { return P.multibody().body_begin; })
        .def_property_readonly("T", &LinearScatteringProblem::T, py::return_value_policy::reference_internal)
        .def_property_readonly("unknowns", &LinearScatteringProblem::unknowns)
        .def_property_readonly("hmatrix_bytes", &LinearScatteringProblem::hmatrix_bytes)
        .def_property_readonly("near_pairs", &LinearScatteringProblem::near_pairs);
}

}  // namespace cbem::py_bind
