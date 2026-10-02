// Python-Anbindung: gekruemmte (quadratische) Elemente mit unstetig linearen Dichten (v0.47, Stufe 2b).
// Spuren: 24 Koeffizienten je Element wie bei den linearen Dichten (h.reshape(-1, 3, 8)).
#include "common.hpp"

#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/curved_problem.hpp"

namespace cbem::py_bind {

namespace {

std::vector<cplx> to_trace24(const CArr& a, std::size_t N, const char* what) {
    if (static_cast<std::size_t>(a.size()) != 24 * N)
        throw py::value_error(std::string(what) + ": erwartet 24 N = " + std::to_string(24 * N) + " Koeffizienten, erhalten " + std::to_string(a.size()));
    return to_cvec(a);
}

py::array_t<double> mids_to_numpy(const QuadraticMesh& m) {
    py::array_t<double> a({static_cast<py::ssize_t>(m.size()), py::ssize_t(3), py::ssize_t(3)});
    double* d = a.mutable_data();
    for (std::size_t t = 0; t < m.size(); ++t)
        for (int e = 0; e < 3; ++e) { d[9 * t + 3 * e] = m.mid[t][e].x; d[9 * t + 3 * e + 1] = m.mid[t][e].y; d[9 * t + 3 * e + 2] = m.mid[t][e].z; }
    return a;
}

}  // namespace

void init_curved(py::module_& m) {
    py::class_<QuadraticMesh>(m, "QuadraticMesh", R"doc(
Quadratische (gekruemmte) Dreieckselemente wie Gmsh ElementOrder 2: die Ecken des ebenen Netzes (flat) und je Element die
Kantenmitten M_01, M_12, M_20 (midpoints, Form (N, 3, 3)) auf der Flaeche. Erzeugt mit quadratic_icosphere,
make_quadratic_sphere oder make_quadratic(mesh, project) mit einer Funktion, die die Sehnenmitte auf die Flaeche abbildet.
)doc")
        .def(py::init([](const TriangleMesh& flat, const RArr& mids) {
                 if (mids.ndim() != 3 || static_cast<std::size_t>(mids.shape(0)) != flat.size() || mids.shape(1) != 3 || mids.shape(2) != 3)
                     throw py::value_error("midpoints: Form (N, 3, 3) erwartet (je Element M_01, M_12, M_20)");
                 QuadraticMesh q; q.flat = flat; q.mid.resize(flat.size());
                 const double* d = mids.data();
                 for (std::size_t t = 0; t < flat.size(); ++t)
                     for (int e = 0; e < 3; ++e) q.mid[t][e] = Vec3(d[9 * t + 3 * e], d[9 * t + 3 * e + 1], d[9 * t + 3 * e + 2]);
                 q.compute_geometry();
                 return q;
             }),
             "flat"_a, "midpoints"_a, "aus ebenem Netz und Kantenmitten; benachbarte Elemente muessen dieselben Mitten verwenden")
        .def_property_readonly("flat", [](const QuadraticMesh& q) { return q.flat; }, "ebenes Netz (Ecken, Konnektivitaet), Kopie")
        .def_property_readonly("midpoints", &mids_to_numpy)
        .def_property_readonly("areas", [](const QuadraticMesh& q) { return to_numpy(q.area); }, "Flaechen der gekruemmten Elemente")
        .def_property_readonly("bulge", [](const QuadraticMesh& q) { return to_numpy(q.bulge); }, "Abstand der Kantenmitten von den Sehnenmitten")
        .def("__len__", &QuadraticMesh::size)
        .def_property_readonly("unknowns", [](const QuadraticMesh& q) { return 24 * q.size(); })
        .def("volume", [](const QuadraticMesh& q, int sub) { return signed_volume(q, sub); }, "sub"_a = 3)
        .def("point", [](const QuadraticMesh& q, std::size_t t, const Vec3& lam) {
                 if (t >= q.size()) throw py::index_error("Elementindex");
                 return q.X(t, {lam.x, lam.y, lam.z});
             }, "t"_a, "lam"_a, "Punkt X(lambda) des Elements t (baryzentrische Parameter)")
        .def("normal", [](const QuadraticMesh& q, std::size_t t, const Vec3& lam) {
                 if (t >= q.size()) throw py::index_error("Elementindex");
                 Vec3 n; q.jacobian(t, {lam.x, lam.y, lam.z}, &n); return n;
             }, "t"_a, "lam"_a)
        .def("__repr__", [](const QuadraticMesh& q) { return "QuadraticMesh(" + std::to_string(q.size()) + " quadratische Elemente)"; });
    m.def("quadratic_icosphere", &quadratic_icosphere, "n"_a, "radius"_a = 1.0, "Ikosaederkugel mit quadratischen Elementen (20 n^2)");
    m.def("make_quadratic_sphere", &make_quadratic_sphere, "mesh"_a, "center"_a, "radius"_a, "Kantenmitten radial auf die Kugel");
    m.def("make_quadratic", [](const TriangleMesh& flat, py::function project) {
              return make_quadratic(flat, [&](const Vec3& p) { return project(p).cast<Vec3>(); });
          }, "mesh"_a, "project"_a, "Kantenmitten = project(Sehnenmitte), project: Funktion (3,) -> (3,)");
    m.def("translated", [](const QuadraticMesh& q, const Vec3& s, real f) { return translated(q, s, f); }, "mesh"_a, "shift"_a, "scale"_a = 1.0);
    m.def("read_gmsh_quadratic", [](const std::string& path, real scale) {
              auto bodies = read_gmsh_quadratic(path, scale);
              py::list out;
              for (auto& b : bodies) out.append(py::make_tuple(b.tag, std::move(b.mesh)));
              return out;
          }, "path"_a, "scale"_a = 1.0,
          "Gmsh 2.2/4.1 zweiter Ordnung (6-Knoten-Dreiecke, Mesh.ElementOrder = 2): Liste von (tag, QuadraticMesh), nach aussen orientiert");
    m.def("write_gmsh22_quadratic", &write_gmsh22_quadratic, "path"_a, "bodies"_a, "Gmsh 2.2 zweiter Ordnung, physikalische Gruppen 1, 2, ...");

    ParamClass<CurvedNearParams>(m, "CurvedNearParams", "Nahquadratur gekruemmter Elemente (doppelt adaptiv)")
        .field("outer_ratio", &CurvedNearParams::outer_ratio, "aeusseres Teilstueck: Umkreisradius < outer_ratio * Abstand")
        .field("inner_ratio", &CurvedNearParams::inner_ratio, "inneres Teilstueck: Umkreisradius < inner_ratio * Abstand")
        .field("outer_depth", &CurvedNearParams::outer_depth, "maximale Halbierungstiefe aussen")
        .field("inner_depth", &CurvedNearParams::inner_depth, "maximale Halbierungstiefe innen");

    py::class_<CurvedKernelEntries>(m, "CurvedKernelEntries", "Eintraege gekruemmter Elemente: je Elementpaar (3, 3, 7), Komponenten [1, e1, e2, e3, e12, e13, e23]")
        .def(py::init<const QuadraticMesh&, cplx, EntryParams, CurvedNearParams>(), "mesh"_a, "k"_a, "params"_a = EntryParams{},
             "near"_a = CurvedNearParams{}, py::keep_alive<1, 2>(), py::call_guard<py::gil_scoped_release>())
        .def("block", [](const CurvedKernelEntries& E, std::size_t i, std::size_t j) {
                 const CurvedBlock B = E.block(i, j);
                 py::array_t<cplx> a({py::ssize_t(3), py::ssize_t(3), py::ssize_t(kCurvedComps)});
                 for (int p = 0; p < 9; ++p) for (int c = 0; c < kCurvedComps; ++c) a.mutable_data()[kCurvedComps * p + c] = B[p][c];
                 return a;
             }, "i"_a, "j"_a)
        .def_property_readonly("near_pairs", &CurvedKernelEntries::near_pairs)
        .def_property_readonly("near_seconds", &CurvedKernelEntries::near_seconds);
    py::class_<CurvedHMatrix>(m, "CurvedHMatrix")
        .def(py::init<const CurvedKernelEntries&, HMatrixParams>(), "entries"_a, "params"_a = HMatrixParams{}, py::keep_alive<1, 2>(),
             py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("stats", &CurvedHMatrix::stats);
    py::class_<CurvedCauchyOperator, BoundaryOperator>(m, "CurvedCauchyOperator")
        .def(py::init<const QuadraticMesh&, const CurvedHMatrix&>(), "mesh"_a, "H"_a, py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

    m.def("plane_wave_trace_curved", [](const QuadraticMesh& q, const Medium& md, real omega, const Vec3& d, const CVec3& p, int sub) {
              if (std::abs(md.chi) > 0) throw py::value_error("gekruemmte Elemente: chirales Aussenmedium noch nicht unterstuetzt");
              return to_numpy(nogil([&] { return project_plane_wave_curved(q, md.k(omega), md.eps, d, p, sub); }));
          }, "mesh"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "sub"_a = 2);
    m.def("extinction_cross_section_curved", [](const QuadraticMesh& q, const CArr& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
              auto h = to_trace24(hs, q.size(), "h_scat");
              return nogil([&] { return extinction_cross_section_curved(q, h, k, eps, d, p); });
          }, "mesh"_a, "h_scat"_a, "k"_a, "eps"_a, "d"_a, "p"_a);
    m.def("far_field_curved", [](const QuadraticMesh& q, const CArr& hs, cplx k, const Vec3& xhat, int sub) {
              auto h = to_trace24(hs, q.size(), "h_scat");
              return nogil([&] { return far_field_curved(q, h, k, xhat, sub); });
          }, "mesh"_a, "h_scat"_a, "k"_a, "xhat"_a, "sub"_a = 2);

    py::class_<CurvedScatteringProblem>(m, "CurvedScatteringProblem", R"doc(
Streuproblem auf quadratischen (gekruemmten) Elementen mit unstetig linearen Dichten (24 Unbekannte je Element): ein oder
mehrere Koerper (auch chiral), achirales Aussenmedium. Kugel mit 320 Elementen: Fehler gegen Mie etwa 0,01 % (eben,
konstant: 5,7 % bei Glas, 7,3 % bei Gold); Konvergenz etwa O(h^4) (docs/results_curved.md).
)doc")
        .def(py::init([](const std::vector<QuadraticMesh>& bodies, const std::vector<Medium>& media, real omega, const Medium& outer,
                         const HMatrixParams& hp, const EntryParams& ep, const CurvedNearParams& np) {
                 if (bodies.empty() || bodies.size() != media.size()) throw py::value_error("je Koerper ein Medium erwartet");
                 return nogil([&] { return std::make_unique<CurvedScatteringProblem>(bodies, media, omega, outer, hp, ep, np); });
             }),
             "bodies"_a, "media"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = HMatrixParams{}, "entries"_a = EntryParams{},
             "near"_a = CurvedNearParams{})
        .def(py::init([](const QuadraticMesh& body, const Medium& medium, real omega, const Medium& outer, const HMatrixParams& hp,
                         const EntryParams& ep, const CurvedNearParams& np) {
                 return nogil([&] {
                     return std::make_unique<CurvedScatteringProblem>(std::vector<QuadraticMesh>{body}, std::vector<Medium>{medium}, omega, outer, hp, ep, np);
                 });
             }),
             "body"_a, "medium"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = HMatrixParams{}, "entries"_a = EntryParams{},
             "near"_a = CurvedNearParams{})
        .def("solve_plane_wave", &CurvedScatteringProblem::solve_plane_wave, "d"_a, "p"_a, "options"_a = SolveOptions{},
             py::call_guard<py::gil_scoped_release>())
        .def("solve_rhs", [](const CurvedScatteringProblem& P, const CArr& b, const SolveOptions& o) {
                 auto bv = to_trace24(b, P.mesh().size(), "b");
                 return nogil([&] { return P.solve_rhs(bv, o); });
             }, "b"_a, "options"_a = SolveOptions{})
        .def_property_readonly("mesh", &CurvedScatteringProblem::mesh, py::return_value_policy::reference_internal)
        .def_property_readonly("body_begin", &CurvedScatteringProblem::body_begin)
        .def_property_readonly("unknowns", &CurvedScatteringProblem::unknowns)
        .def_property_readonly("hmatrix_bytes", &CurvedScatteringProblem::hmatrix_bytes)
        .def_property_readonly("near_seconds", &CurvedScatteringProblem::near_seconds);
}

}  // namespace cbem::py_bind
