// Python-Anbindung: gekruemmte (quadratische) Elemente mit unstetig linearen Dichten (v0.47, Stufe 2b).
// Spuren: 24 Koeffizienten je Element wie bei den linearen Dichten (h.reshape(-1, 3, 8)).
#include "common.hpp"

#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/curved_layered_problem.hpp"
#include "cbem/problems/curved_problem.hpp"
#include "cbem/sources/chiral_incidence.hpp"
#include "cbem/sources/curved_near_field.hpp"

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

template <class F> py::array_t<cplx> apply_checked(std::size_t n, const CArr& x, F&& f) {
    if (static_cast<std::size_t>(x.size()) != n)
        throw py::value_error("x: Laenge " + std::to_string(n) + " erwartet, erhalten " + std::to_string(x.size()));
    auto xv = to_cvec(x);
    std::vector<cplx> y;
    nogil([&] { f(xv, y); });
    return to_numpy(y);
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

    ParamClass<CurvedNearParams>(m, "CurvedNearParams", R"doc(
Nahquadratur gekruemmter Elemente fuer nahe, getrennte Paare. Voreinstellung (v0.52): Singularitaetssubtraktion am
Tangentialdreieck, Gauss 5 x 5 an den Blaettern (outer_rule = correction_rule = 5), Kriterien 1,5/1,5 (Fehler 0,6-2,6e-6,
3-mal schneller als v0.49). Schnell: Regeln 4, Kriterien 1.0 (Fehler 2,5-4e-5). outer_rule/correction_rule = 0: Dunavant 7
(v0.49 mit 0,3/0,3). subtract = False: doppelt adaptiv (outer_ratio, inner_ratio) als Referenz.
)doc")
        .field("subtract", &CurvedNearParams::subtract, "Singularitaetssubtraktion am Tangentialdreieck (sonst doppelt adaptiv)")
        .field("subtract_outer_ratio", &CurvedNearParams::subtract_outer_ratio, "mit Subtraktion: aeusseres Kriterium")
        .field("correction_ratio", &CurvedNearParams::correction_ratio, "mit Subtraktion: Kriterium fuer den Rest")
        .field("outer_rule", &CurvedNearParams::outer_rule, "mit Subtraktion: Blattregel aussen (0 = Dunavant 7, n = Gauss n x n)")
        .field("correction_rule", &CurvedNearParams::correction_rule, "mit Subtraktion: Blattregel der Korrektur (0 = Dunavant 7, n = Gauss n x n)")
        .field("ss_orders", &CurvedNearParams::ss_orders,
               "Sauter-Schwab [Ecke, Kante, Selbstterm] je Richtung (xi, eta1, eta2, eta3); erste Zahl 0: isotrop mit EntryParams.ss_order")
        .field("outer_ratio", &CurvedNearParams::outer_ratio, "doppelt adaptiv: aeusseres Teilstueck, Umkreisradius < outer_ratio * Abstand")
        .field("inner_ratio", &CurvedNearParams::inner_ratio, "doppelt adaptiv: inneres Teilstueck, Umkreisradius < inner_ratio * Abstand")
        .field("outer_depth", &CurvedNearParams::outer_depth, "maximale Halbierungstiefe aussen")
        .field("inner_depth", &CurvedNearParams::inner_depth, "maximale Halbierungstiefe innen")
        .field("adapt_to_boundary", &CurvedNearParams::adapt_to_boundary,
               "duenne Schichten (v0.62): Randabstand fuer Teilstuecke auf einer Seite der Flaeche, Abstaende zur gekruemmten Flaeche");

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
        .def(py::init<const CurvedKernelEntries&, HMatrixParams>(), "entries"_a, "params"_a = curved_hmatrix_params(), py::keep_alive<1, 2>(),
             py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("stats", &CurvedHMatrix::stats);
    py::class_<CurvedCauchyOperator, BoundaryOperator>(m, "CurvedCauchyOperator")
        .def(py::init<const QuadraticMesh&, const CurvedHMatrix&>(), "mesh"_a, "H"_a, py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

    m.def("plane_wave_trace_curved", [](const QuadraticMesh& q, const Medium& md, real omega, const Vec3& d, const CVec3& p, int sub) {
              const cplx k = plane_wave_incidence(md, omega, d, p).k;     // chiral: Helizitaetswelle mit k_sigma (v0.59)
              return to_numpy(nogil([&] { return project_plane_wave_curved(q, k, md.eps, d, p, sub); }));
          }, "mesh"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "sub"_a = 2);
    m.def("extinction_cross_section_curved", [](const QuadraticMesh& q, const CArr& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
              auto h = to_trace24(hs, q.size(), "h_scat");
              return nogil([&] { return extinction_cross_section_curved(q, h, k, eps, d, p); });
          }, "mesh"_a, "h_scat"_a, "k"_a, "eps"_a, "d"_a, "p"_a);
    m.def("far_field_curved", [](const QuadraticMesh& q, const CArr& hs, cplx k, const Vec3& xhat, int sub) {
              auto h = to_trace24(hs, q.size(), "h_scat");
              return nogil([&] { return far_field_curved(q, h, k, xhat, sub); });
          }, "mesh"_a, "h_scat"_a, "k"_a, "xhat"_a, "sub"_a = 2);

    // Nahfeld und Kraefte (v0.58)
    m.def("project_incident_curved", [](const QuadraticMesh& q, const IncidentField& inc, const Medium& md, int sub) {
              if (sub < 1) throw py::value_error("sub >= 1");
              return to_numpy(nogil([&] { return project_incident_curved(q, inc, md, sub); }));
          }, "mesh"_a, "incident"_a, "medium"_a, "sub"_a = 2, "Spur eines einfallenden Feldes in der psi-Basis (24 N)");
    m.def("scattered_field_curved", [](const QuadraticMesh& q, const CArr& hs, cplx k, const RArr& pts) {
              auto h = to_trace24(hs, q.size(), "h_scat");
              auto x = to_points(pts);
              return multivectors_to_numpy(nogil([&] { return scattered_field_curved(q, h, k, x); }));
          }, "mesh"_a, "h_scat"_a, "k"_a, "points"_a, "Streufeld-Multivektoren (M, 8) aus der Streuspur (24 N); direkte Summation");
    m.def("scattered_field_curved_hmatrix", [](const QuadraticMesh& q, const CArr& hs, cplx k, const RArr& pts, real eps) {
              if (!(eps > 0)) throw py::value_error("eps > 0");
              auto h = to_trace24(hs, q.size(), "h_scat");
              auto x = to_points(pts);
              return multivectors_to_numpy(nogil([&] { return scattered_field_curved_hmatrix(q, h, k, x, eps); }));
          }, "mesh"_a, "h_scat"_a, "k"_a, "points"_a, "eps"_a = 1e-6, R"doc(
Streufeld wie scattered_field_curved mit einer H-Matrix (v0.61; ACA+ auf dem Dirac-Kern an den Quadraturpunkten, Bloecke sofort
angewandt): 5120 Elemente x 40 000 Punkte 3,7 s statt 27 s, Fehler etwa 3e-7 bezogen auf max |F_s| bei eps = 1e-6.
)doc");
    m.def("curved_near_field_options", &curved_near_field_options,
          "Voreinstellung der Nahfeldauswertung gekruemmter Elemente: H-Matrix ab 4000 Punkten, eps = 1e-6 (eben 2000 Punkte, 1e-4)");
    m.def("exterior_near_field_curved", [](const QuadraticMesh& q, const CArr& h, const Medium& md, real omega, const Vec3& d, const CVec3& p,
                                           const RArr& pts, const NearFieldOptions& o) {
              auto hv = to_trace24(h, q.size(), "h");
              auto x = to_points(pts);
              return near_field_to_dict(nogil([&] { return exterior_near_field_curved(q, hv, md, omega, d, p, x, o); }));
          }, "outer"_a, "h"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "points"_a, "options"_a = curved_near_field_options(), R"doc(
Gesamtes Nahfeld auf gekruemmten Elementen bei Anregung durch eine ebene Welle; h = Gesamtspur (24 N). dict wie
exterior_near_field (E, H, inside, too_close, enhancement, chirality). Goldkugel 1280 Elemente: |E|^2 auf 0,05 % (eben 2-6 %).
Ab options.hmatrix_min_points Punkten (4000) mit H-Matrix (ACA-Toleranz options.eps, Voreinstellung 1e-6; v0.61).
)doc");
    m.def("exterior_near_field_curved", [](const QuadraticMesh& q, const CArr& h, const CArr& b, const Medium& md, real omega,
                                           const IncidentField& inc, const RArr& pts, const NearFieldOptions& o) {
              auto hv = to_trace24(h, q.size(), "h");
              auto bv = to_trace24(b, q.size(), "b");
              auto x = to_points(pts);
              return near_field_to_dict(nogil([&] { return exterior_near_field_curved(q, hv, bv, md, omega, inc, x, o); }));
          }, "outer"_a, "h"_a, "b"_a, "medium"_a, "omega"_a, "incident"_a, "points"_a, "options"_a = curved_near_field_options(),
          "allgemeine Anregung: b = project_incident_curved(outer, incident, medium)");
    m.def("near_field_evaluator_curved", [](const QuadraticMesh& q, const CArr& h, const CArr& b, const Medium& md, real omega,
                                            std::shared_ptr<IncidentField> inc, const NearFieldOptions& o) {
              struct Data { QuadraticMesh q; std::vector<cplx> h; };
              auto data = std::make_shared<const Data>(Data{q, to_trace24(h, q.size(), "h")});
              PyNearFieldEval e;
              e.owner = data;
              e.f = make_near_field_eval_curved(data->q, data->h, to_trace24(b, q.size(), "b"), md, omega, inc, o);
              return e;
          }, "outer"_a, "h"_a, "b"_a, "medium"_a, "omega"_a, "incident"_a, "options"_a = curved_near_field_options(),
          "Feldauswerter fuer Kraefte und Gradienten (force_on_sphere, force_on_offset, fields_with_gradients) auf gekruemmten Elementen");

    // einfallende Felder aus Python (v0.60): Projektion in zwei Schritten (CustomField.project waehlt sie fuer QuadraticMesh),
    // Nahfeld und Feldauswerter wie im ebenen Pfad (Streufeld ohne GIL, fields(x) ein vektorisierter Aufruf je Auswertung)
    m.def("_quadrature_points_curved", [](const QuadraticMesh& q, int sub) {
              if (sub < 1) throw py::value_error("sub >= 1");
              return points_to_numpy(projection_points_curved(q, sub));
          }, "mesh"_a, "sub"_a = 2, "Quadraturpunkte der Projektion auf gekruemmte Elemente (N q, 3), elementweise");
    m.def("_project_samples_curved", [](const QuadraticMesh& q, const CArr& Ea, const CArr& Ha, const Medium& md, int sub) {
              if (sub < 1) throw py::value_error("sub >= 1");
              const std::size_t M = projection_points_curved(q, sub).size();
              auto take = [&](const CArr& a, const char* what) {
                  if (a.ndim() != 2 || static_cast<std::size_t>(a.shape(0)) != M || a.shape(1) != 3)
                      throw py::value_error(std::string(what) + ": Form (" + std::to_string(M) + ", 3) erwartet");
                  std::vector<CVec3> v(M); const cplx* d = a.data();
                  for (std::size_t i = 0; i < M; ++i) v[i] = CVec3{d[3 * i], d[3 * i + 1], d[3 * i + 2]};
                  return v;
              };
              auto E = take(Ea, "E"), H = take(Ha, "H");
              return to_numpy(nogil([&] { return project_samples_curved(q, E, H, md, sub); }));
          }, "mesh"_a, "E"_a, "H"_a, "medium"_a, "sub"_a = 2, "Spur (24 N) aus Feldwerten an den Quadraturpunkten");
    m.def("exterior_near_field_curved", [](const QuadraticMesh& q, const CArr& h, const CArr& b, const Medium& md, real omega,
                                           py::object incident, const RArr& pts, const NearFieldOptions& o) {
              check_python_field(incident);
              auto hv = to_trace24(h, q.size(), "h");
              auto bv = to_trace24(b, q.size(), "b");
              auto x = to_points(pts);
              return near_field_to_dict(nogil([&] {
                  return add_python_incident(exterior_near_field_curved(q, hv, bv, md, omega, ZeroField(), x, o), incident, x);
              }));
          }, "outer"_a, "h"_a, "b"_a, "medium"_a, "omega"_a, "incident"_a, "points"_a, "options"_a = curved_near_field_options(),
          "einfallendes Feld aus Python (Objekt mit fields(x) -> (E, H), z. B. CustomField); b = incident.project(outer, medium)");
    m.def("near_field_evaluator_curved", [](const QuadraticMesh& q, const CArr& h, const CArr& b, const Medium& md, real omega,
                                            py::object incident, const NearFieldOptions& o) {
              check_python_field(incident);
              struct Data { QuadraticMesh q; std::vector<cplx> h, b; };
              auto data = std::make_shared<const Data>(Data{q, to_trace24(h, q.size(), "h"), to_trace24(b, q.size(), "b")});
              auto field = hold_python_object(incident);
              PyNearFieldEval e;
              e.owner = data;
              // wird von den Kraftfunktionen ohne GIL aufgerufen, je Auswertung einmal mit allen Punkten
              e.f = [data, md, omega, field, o](const std::vector<Vec3>& pts) {
                  return add_python_incident(exterior_near_field_curved(data->q, data->h, data->b, md, omega, ZeroField(), pts, o), *field, pts);
              };
              return e;
          }, "outer"_a, "h"_a, "b"_a, "medium"_a, "omega"_a, "incident"_a, "options"_a = curved_near_field_options(),
          "Feldauswerter fuer ein einfallendes Feld aus Python auf gekruemmten Elementen");

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
             "bodies"_a, "media"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = curved_hmatrix_params(), "entries"_a = EntryParams{},
             "near"_a = CurvedNearParams{})
        .def(py::init([](const QuadraticMesh& body, const Medium& medium, real omega, const Medium& outer, const HMatrixParams& hp,
                         const EntryParams& ep, const CurvedNearParams& np) {
                 return nogil([&] {
                     return std::make_unique<CurvedScatteringProblem>(std::vector<QuadraticMesh>{body}, std::vector<Medium>{medium}, omega, outer, hp, ep, np);
                 });
             }),
             "body"_a, "medium"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = curved_hmatrix_params(), "entries"_a = EntryParams{},
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

    // geschichtete Koerper (v0.62)
    m.def("offset_surface", [](const QuadraticMesh& q, real d) { return nogil([&] { return offset_surface(q, d); }); }, "mesh"_a, "d"_a,
          "Parallelflaeche eines quadratischen Netzes (glatte Flaechen; Ecken und Kantenmitten entlang der gemittelten Normalen)");
    m.def("curved_layered_near_params", &curved_layered_near_params,
          "Nahquadratur fuer geschichtete Koerper auf gekruemmten Elementen (adapt_to_boundary = True)");
    py::class_<CurvedLayeredGeometry>(m, "CurvedLayeredGeometry", R"doc(
Grenzflaechengraph aus quadratischen Flaechen (wie LayeredGeometry): geschlossene Flaechen, je mit Innen- und Aussengebiet;
Gebiet 0 ist der Aussenraum.
)doc")
        .def(py::init<Medium>(), "exterior"_a = Medium{})
        .def("add_region", &CurvedLayeredGeometry::add_region, "medium"_a, "neues Gebiet; gibt seinen Index zurueck")
        .def("add_surface", &CurvedLayeredGeometry::add_surface, "mesh"_a, "inside"_a, "outside"_a, "Flaeche zwischen zwei Gebieten")
        .def("add_body", [](CurvedLayeredGeometry& g, const QuadraticMesh& s, const Medium& in, int parent) { return add_body(g, s, in, parent); },
             "surface"_a, "inner"_a, "parent"_a = 0, "homogener Koerper; gibt sein Gebiet zurueck")
        .def("add_layered_body", [](CurvedLayeredGeometry& g, const std::vector<QuadraticMesh>& s, const std::vector<Medium>& media, int parent) {
                 if (s.size() != media.size()) throw py::value_error("je Flaeche ein Medium erwartet");
                 return add_layered_body(g, s, media, parent);
             }, "surfaces_outer_first"_a, "media"_a, "parent"_a = 0, "verschachtelte Flaechen von aussen nach innen")
        .def("add_coated_body", [](CurvedLayeredGeometry& g, const QuadraticMesh& s, const Medium& core, const std::vector<Coating>& c, bool outward,
                                   int parent) { return nogil([&] { return add_coated_body(g, s, core, c, outward, parent); }); },
             "surface"_a, "core"_a, "coatings"_a, "outward"_a = true, "parent"_a = 0,
             "beschichteter Koerper, Schichten von innen nach aussen (Parallelflaechen, offset_surface)")
        // Vektorfelder als Kopien (Referenzen in die Vektoren wuerden nach add_surface ungueltig)
        .def_property_readonly("surfaces", [](const CurvedLayeredGeometry& g) { return g.surfaces; })
        .def_property_readonly("inside", [](const CurvedLayeredGeometry& g) { return g.inside; })
        .def_property_readonly("outside", [](const CurvedLayeredGeometry& g) { return g.outside; })
        .def_property_readonly("region_medium", [](const CurvedLayeredGeometry& g) { return g.region_medium; })
        .def_property_readonly("elements", &CurvedLayeredGeometry::elements);

    py::class_<CurvedLayeredTransmissionOperator>(m, "CurvedLayeredTransmissionOperator")
        .def("apply", [](const CurvedLayeredTransmissionOperator& T, const CArr& x) {
            return apply_checked(T.size(), x, [&](const std::vector<cplx>& a, std::vector<cplx>& y) { T.apply(a, y); });
        }, "x"_a)
        .def("precondition", [](const CurvedLayeredTransmissionOperator& T, const CArr& x) {
            return apply_checked(T.size(), x, [&](const std::vector<cplx>& a, std::vector<cplx>& y) { T.precondition(a, y); });
        }, "x"_a)
        .def("__len__", &CurvedLayeredTransmissionOperator::size);

    py::class_<CurvedLayeredScatteringProblem>(m, "CurvedLayeredScatteringProblem", R"doc(
Geschichtete und beschichtete Koerper auf gekruemmten Elementen (v0.62). Goldkern mit Glasschale, 2 x 1280 Elemente: sigma_ext
gegen Aden-Kerker 2e-6 (d = 0,2) bis 2,5e-5 (duenne Schichten bis d = 0,01), die Wirkung der Schicht direkt gegen die Rechnung
ohne Schicht auf 2e-4 (konstante Dichten 0,3-0,9 % bzw. 1-2 % gegen eine neutrale Vergleichsrechnung). Duenne Schichten
kosten Aufbauzeit (d = 0,01: etwa 440 s mit 12 Threads).
)doc")
        .def(py::init([](const CurvedLayeredGeometry& g, real omega, const HMatrixParams& hp, const EntryParams& ep, const CurvedNearParams& np) {
                 return nogil([&] { return std::make_unique<CurvedLayeredScatteringProblem>(g, omega, hp, ep, np); });
             }),
             "geometry"_a, "omega"_a, "hmatrix"_a = curved_hmatrix_params(), "entries"_a = EntryParams{}, "near"_a = curved_layered_near_params())
        .def("solve_plane_wave", &CurvedLayeredScatteringProblem::solve_plane_wave, "d"_a, "p"_a, "options"_a = SolveOptions{},
             py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("mesh", &CurvedLayeredScatteringProblem::mesh, py::return_value_policy::reference_internal, "alle Flaechen")
        .def("surface_begin", &CurvedLayeredScatteringProblem::surface_begin, "s"_a)
        .def_property_readonly("T", &CurvedLayeredScatteringProblem::T, py::return_value_policy::reference_internal)
        .def_property_readonly("unknowns", &CurvedLayeredScatteringProblem::unknowns)
        .def_property_readonly("hmatrix_bytes", &CurvedLayeredScatteringProblem::hmatrix_bytes)
        .def_property_readonly("near_seconds", &CurvedLayeredScatteringProblem::near_seconds);
}

}  // namespace cbem::py_bind
