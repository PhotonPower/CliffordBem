// Python-Anbindung: Streuprobleme (homogen, geschichtet, Duennschicht-Naeherung, Zweitor) (v0.42).
//
// Aufbau und Loesen laufen ohne GIL; die Problemklassen kopieren ihre Netze, Python-Objekte duerfen danach freigegeben werden.
#include "common.hpp"

#include <sstream>

#include "cbem/problems/layered_problem.hpp"
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/problems/thin_layer_problem.hpp"
#include "cbem/problems/twoport_layer_problem.hpp"

namespace cbem::py_bind {

namespace {

template <class R> std::string result_repr(const char* name, const R& r) {
    std::ostringstream o;
    o.precision(10);
    o << name << "(sigma_ext=" << r.sigma_ext << ", forward=" << py::repr(py::cast(r.forward)).template cast<std::string>()
      << ", iterations=" << r.iterations << ", residual=" << r.residual << ", len(h)=" << r.h.size() << ")";
    return o.str();
}

template <class F> py::array_t<cplx> apply_checked(std::size_t n, const CArr& x, F&& f) {
    if (static_cast<std::size_t>(x.size()) != n)
        throw py::value_error("x: Laenge " + std::to_string(n) + " erwartet, erhalten " + std::to_string(x.size()));
    auto xv = to_cvec(x);
    std::vector<cplx> y;
    nogil([&] { f(xv, y); });
    return to_numpy(y);
}

void check_media(std::size_t bodies, std::size_t media) {
    if (bodies == 0) throw py::value_error("mindestens ein Koerper erwartet");
    if (bodies != media) throw py::value_error("je Koerper ein Medium erwartet (" + std::to_string(bodies) + " Koerper, " +
                                               std::to_string(media) + " Medien)");
}

}  // namespace

void init_problems(py::module_& m) {
    // --- Ergebnisse -------------------------------------------------------------------------------------------------------------
    py::class_<PlaneWaveResult>(m, "PlaneWaveResult", "Loesung von ScatteringProblem")
        .def_readonly("sigma_ext", &PlaneWaveResult::sigma_ext, "Extinktionsquerschnitt (optisches Theorem), in Einheiten^2")
        .def_readonly("forward", &PlaneWaveResult::forward, "Vorwaertsamplitude S(0) (Bohren-Huffman)")
        .def_readonly("iterations", &PlaneWaveResult::iterations)
        .def_readonly("residual", &PlaneWaveResult::residual)
        .def_property_readonly("h", [](const PlaneWaveResult& r) { return to_numpy(r.h); }, "Gesamtspur aussen (8 N), Kopie")
        .def("__repr__", [](const PlaneWaveResult& r) { return result_repr("PlaneWaveResult", r); });
    py::class_<LayeredResult>(m, "LayeredResult", "Loesung der geschichteten Probleme")
        .def_readonly("sigma_ext", &LayeredResult::sigma_ext)
        .def_readonly("forward", &LayeredResult::forward)
        .def_readonly("iterations", &LayeredResult::iterations)
        .def_readonly("residual", &LayeredResult::residual)
        .def_property_readonly("h", [](const LayeredResult& r) { return to_numpy(r.h); },
                               "Spuren (Layered: alle Flaechen; Zweitor: Gesamtvektor, Aussenspur = erste 8 N), Kopie")
        .def("__repr__", [](const LayeredResult& r) { return result_repr("LayeredResult", r); });

    py::class_<Coating>(m, "Coating", "Schicht mit Dicke und Medium")
        .def(py::init([](real d, const Medium& me) { return Coating{d, me}; }), "thickness"_a, "medium"_a)
        .def_readwrite("thickness", &Coating::thickness)
        .def_readwrite("medium", &Coating::medium)
        .def("__repr__", [](const Coating& c) {
            return "Coating(thickness=" + py::repr(py::float_(c.thickness)).cast<std::string>() + ", medium=" +
                   py::repr(py::cast(c.medium)).cast<std::string>() + ")";
        });

    // --- ScatteringProblem ------------------------------------------------------------------------------------------------------
    py::class_<ScatteringProblem>(m, "ScatteringProblem", R"doc(
Streuproblem fuer einen oder mehrere homogene Koerper (auch chiral) in einem Aussenmedium (auch chiral).
Baut die Randoperatoren als H-Matrizen und T_1; loest fuer ebene Wellen oder beliebige rechte Seiten.

Laengen in einer frei gewaehlten Einheit L, omega = k0 = 2 pi L / lambda (eps0 = mu0 = 1); sigma_ext in L^2.

Beispiel::

    P = ScatteringProblem(make_icosphere(8), Medium(eps=-11 + 1.2j), omega=0.5)
    r = P.solve_plane_wave(d=(0, 0, 1), p=(1, 0, 0))
    Q = r.sigma_ext / pi
)doc")
        .def(py::init([](const std::vector<TriangleMesh>& bodies, const std::vector<Medium>& media, real omega, const Medium& outer,
                         const HMatrixParams& hp, const EntryParams& ep, bool union_interior) {
                 check_media(bodies.size(), media.size());
                 return nogil([&] { return std::make_unique<ScatteringProblem>(bodies, media, omega, outer, hp, ep, union_interior); });
             }),
             "bodies"_a, "media"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = HMatrixParams{}, "entries"_a = EntryParams{},
             "union_interior"_a = false, "mehrere Koerper (getrennte, nach aussen orientierte geschlossene Flaechen), Medium je Koerper")
        .def(py::init([](const TriangleMesh& body, const Medium& medium, real omega, const Medium& outer, const HMatrixParams& hp,
                         const EntryParams& ep) {
                 return nogil([&] {
                     return std::make_unique<ScatteringProblem>(std::vector<TriangleMesh>{body}, std::vector<Medium>{medium}, omega, outer, hp, ep);
                 });
             }),
             "body"_a, "medium"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = HMatrixParams{}, "entries"_a = EntryParams{}, "ein Koerper")
        .def("solve_plane_wave", &ScatteringProblem::solve_plane_wave, "d"_a, "p"_a, "options"_a = SolveOptions{},
             py::call_guard<py::gil_scoped_release>(),
             "ebene Welle mit Richtung d und (komplexer) Polarisation p; im chiralen Aussenmedium nur zirkular")
        .def("solve_rhs", [](const ScatteringProblem& P, const CArr& b, const SolveOptions& o) {
                 auto bv = to_trace(b, P.mesh().size(), "b");
                 return nogil([&] { return P.solve_rhs(bv, o); });
             }, "b"_a, "options"_a = SolveOptions{}, "beliebige einfallende Spur b (z. B. project_dipole); r.h = Gesamtspur")
        .def("apply", [](const ScatteringProblem& P, const CArr& x) {
                 return apply_checked(P.T().size(), x, [&](const std::vector<cplx>& a, std::vector<cplx>& y) { P.T().apply(a, y); });
             }, "x"_a, "T_1 x")
        .def("use_recycling", &ScatteringProblem::use_recycling, "max_recycle"_a = 120, "Krylov-Recycling fuer folgende Loesungen")
        .def("use_gcrodr", &ScatteringProblem::use_gcrodr, "k"_a = 20, "m"_a = 80, "GCRO-DR mit k harmonischen Ritz-Vektoren")
        .def_property_readonly("recycled", &ScatteringProblem::recycled)
        .def("use_block_preconditioner", &ScatteringProblem::use_block_preconditioner, "groups"_a,
             py::call_guard<py::gil_scoped_release>(), "Blockvorkonditionierung (Gruppen aus group_by_clusters / group_by_features)")
        .def("use_hodlr_preconditioner", &ScatteringProblem::use_hodlr_preconditioner, "params"_a = HodlrParams{},
             py::call_guard<py::gil_scoped_release>(), "HODLR-Faktorisierung als Vorkonditionierer")
        .def_property_readonly("hodlr_info", [](const ScatteringProblem& P) -> py::object {
            const HodlrSolver* h = P.hodlr();
            if (!h) return py::none();
            return py::dict("max_rank"_a = h->max_rank(), "bytes"_a = h->bytes(), "seconds"_a = h->seconds());
        })
        .def_property_readonly("block_info", [](const ScatteringProblem& P) -> py::object {
            const BlockPreconditioner* b = P.block_preconditioner();
            if (!b) return py::none();
            return py::dict("max_block"_a = b->max_block(), "fraction"_a = b->fraction(), "seconds"_a = b->seconds());
        })
        .def("system_entry", [](const ScatteringProblem& P, std::size_t i, std::size_t j) { return mat8_to_numpy(P.system_entry(i, j)); },
             "i"_a, "j"_a, "8x8-Block T_ij des Systems (exakte Eintraege)")
        .def("inner_block", [](const ScatteringProblem& P, const std::vector<std::size_t>& B) {
                 Matrix M = nogil([&] { return P.inner_block(B); });
                 return matrix_to_numpy(M);
             }, "triangles"_a, "Innenoperator E_1 auf B x B, dicht (8|B| x 8|B|)")
        .def_property_readonly("mesh", &ScatteringProblem::mesh, py::return_value_policy::reference_internal, "Vereinigung aller Koerper")
        .def_property_readonly("body_begin", [](const ScatteringProblem& P) { return P.multibody().body_begin; },
                               "erster Dreiecksindex je Koerper, zuletzt die Gesamtzahl")
        .def_property_readonly("bodies", [](const ScatteringProblem& P) { return P.multibody().parts; }, "Teilnetze (Kopien)")
        .def_property_readonly("T", &ScatteringProblem::T, py::return_value_policy::reference_internal, "Transmissionsoperator")
        .def_property_readonly("unknowns", [](const ScatteringProblem& P) { return P.T().size(); })
        .def_property_readonly("hmatrix_bytes", &ScatteringProblem::hmatrix_bytes);

    // --- Geschichtete Koerper ---------------------------------------------------------------------------------------------------
    py::class_<LayeredGeometry>(m, "LayeredGeometry", R"doc(
Grenzflaechengraph: geschlossene Flaechen, je mit Innen- und Aussengebiet; Gebiet 0 ist der Aussenraum.
)doc")
        .def(py::init<Medium>(), "exterior"_a = Medium{})
        .def("add_region", &LayeredGeometry::add_region, "medium"_a, "neues Gebiet; gibt seinen Index zurueck")
        .def("add_surface", &LayeredGeometry::add_surface, "mesh"_a, "inside"_a, "outside"_a, "Flaeche zwischen zwei Gebieten")
        .def("add_body", [](LayeredGeometry& g, const TriangleMesh& s, const Medium& in, int parent) { return add_body(g, s, in, parent); },
             "surface"_a, "inner"_a, "parent"_a = 0, "homogener Koerper; gibt sein Gebiet zurueck")
        .def("add_layered_body", [](LayeredGeometry& g, const std::vector<TriangleMesh>& s, const std::vector<Medium>& media, int parent) {
                 if (s.size() != media.size()) throw py::value_error("je Flaeche ein Medium erwartet");
                 return add_layered_body(g, s, media, parent);
             }, "surfaces_outer_first"_a, "media"_a, "parent"_a = 0, "verschachtelte Flaechen von aussen nach innen")
        .def("add_coated_body", [](LayeredGeometry& g, const TriangleMesh& s, const Medium& core, const std::vector<Coating>& c, bool outward,
                                   int parent) { return nogil([&] { return add_coated_body(g, s, core, c, outward, parent); }); },
             "surface"_a, "core"_a, "coatings"_a, "outward"_a = true, "parent"_a = 0,
             "beschichteter Koerper, Schichten von innen nach aussen (Parallelflaechen)")
        // Vektorfelder als Kopien (Referenzen in die Vektoren wuerden nach add_surface ungueltig)
        .def_property_readonly("surfaces", [](const LayeredGeometry& g) { return g.surfaces; })
        .def_property_readonly("inside", [](const LayeredGeometry& g) { return g.inside; })
        .def_property_readonly("outside", [](const LayeredGeometry& g) { return g.outside; })
        .def_property_readonly("region_medium", [](const LayeredGeometry& g) { return g.region_medium; })
        .def_property_readonly("triangles", &LayeredGeometry::triangles);

    py::class_<LayeredTransmissionOperator>(m, "LayeredTransmissionOperator")
        .def("apply", [](const LayeredTransmissionOperator& T, const CArr& x) {
            return apply_checked(T.size(), x, [&](const std::vector<cplx>& a, std::vector<cplx>& y) { T.apply(a, y); });
        }, "x"_a)
        .def("precondition", [](const LayeredTransmissionOperator& T, const CArr& x) {
            return apply_checked(T.size(), x, [&](const std::vector<cplx>& a, std::vector<cplx>& y) { T.precondition(a, y); });
        }, "x"_a)
        .def("__len__", &LayeredTransmissionOperator::size);

    py::class_<LayeredScatteringProblem>(m, "LayeredScatteringProblem", "exakte Rechnung fuer geschichtete/beschichtete Koerper")
        .def(py::init([](const LayeredGeometry& g, real omega, const HMatrixParams& hp, const EntryParams& ep) {
                 return nogil([&] { return std::make_unique<LayeredScatteringProblem>(g, omega, hp, ep); });
             }),
             "geometry"_a, "omega"_a, "hmatrix"_a = HMatrixParams{}, "entries"_a = layered_entry_params())
        .def("solve_plane_wave", &LayeredScatteringProblem::solve_plane_wave, "d"_a, "p"_a, "options"_a = SolveOptions{},
             py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("mesh", &LayeredScatteringProblem::mesh, py::return_value_policy::reference_internal, "alle Flaechen")
        .def("surface_begin", &LayeredScatteringProblem::surface_begin, "s"_a)
        .def_property_readonly("T", &LayeredScatteringProblem::T, py::return_value_policy::reference_internal)
        .def_property_readonly("hmatrix_bytes", &LayeredScatteringProblem::hmatrix_bytes)
        .def_property_readonly("near_pairs", &LayeredScatteringProblem::near_pairs)
        .def_property_readonly("near_seconds", &LayeredScatteringProblem::near_seconds);

    // --- Duennschicht-Naeherung -------------------------------------------------------------------------------------------------
    py::class_<ThinBody>(m, "ThinBody", "Koerper mit duennen Schichten (eine Referenzflaeche)")
        .def(py::init([](const TriangleMesh& s, const Medium& core, const std::vector<Coating>& c, real f) { return ThinBody{s, core, c, f}; }),
             "surface"_a, "core"_a, "coatings"_a = std::vector<Coating>{}, "inner_fraction"_a = 0.0)
        .def_readwrite("surface", &ThinBody::surface)
        .def_readwrite("core", &ThinBody::core)
        .def_property("coatings", [](const ThinBody& b) { return b.coatings; }, [](ThinBody& b, std::vector<Coating> c) { b.coatings = std::move(c); })
        .def_readwrite("inner_fraction", &ThinBody::inner_fraction);

    py::class_<ThinLayerScatteringProblem>(m, "ThinLayerScatteringProblem", "Duennschicht-Naeherung (Dirac-Form 2. Ordnung)")
        .def(py::init([](const TriangleMesh& s, const Medium& core, const std::vector<Coating>& c, real omega, const Medium& outer, real f,
                         const HMatrixParams& hp, const EntryParams& ep, ThinLayerModel model) {
                 return nogil([&] { return std::make_unique<ThinLayerScatteringProblem>(s, core, c, omega, outer, f, hp, ep, model); });
             }),
             "surface"_a, "core"_a, "coatings"_a, "omega"_a, "outer"_a = Medium{}, "inner_fraction"_a = 0.0, "hmatrix"_a = HMatrixParams{},
             "entries"_a = EntryParams{}, "model"_a = ThinLayerModel::Dirac2Fit)
        .def(py::init([](const std::vector<ThinBody>& b, real omega, const Medium& outer, const HMatrixParams& hp, const EntryParams& ep,
                         ThinLayerModel model) {
                 if (b.empty()) throw py::value_error("mindestens ein Koerper erwartet");
                 return nogil([&] { return std::make_unique<ThinLayerScatteringProblem>(b, omega, outer, hp, ep, model); });
             }),
             "bodies"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = HMatrixParams{}, "entries"_a = EntryParams{},
             "model"_a = ThinLayerModel::Dirac2Fit)
        .def("solve_plane_wave", &ThinLayerScatteringProblem::solve_plane_wave, "d"_a, "p"_a, "options"_a = SolveOptions{},
             py::call_guard<py::gil_scoped_release>())
        .def("apply", [](const ThinLayerScatteringProblem& P, const CArr& x) {
            return apply_checked(8 * P.mesh().size(), x, [&](const std::vector<cplx>& a, std::vector<cplx>& y) { P.apply(a, y); });
        }, "x"_a, "T_eff x")
        .def("precondition", [](const ThinLayerScatteringProblem& P, const CArr& x) {
            return apply_checked(8 * P.mesh().size(), x, [&](const std::vector<cplx>& a, std::vector<cplx>& y) { P.precondition(a, y); });
        }, "x"_a)
        .def_property_readonly("mesh", &ThinLayerScatteringProblem::mesh, py::return_value_policy::reference_internal)
        .def("body_begin", &ThinLayerScatteringProblem::body_begin, "b"_a)
        .def_property_readonly("hmatrix_bytes", &ThinLayerScatteringProblem::hmatrix_bytes);

    // --- Zweitor ------------------------------------------------------------------------------------------------------------------
    py::class_<TwoPortBody>(m, "TwoPortBody", "Koerper fuer das Zweitor: Flaechen Gamma_0 (Kern) ... Gamma_L (aussen)")
        .def(py::init([](const std::vector<TriangleMesh>& s, const Medium& core, const std::vector<Coating>& l) {
                 if (s.size() != l.size() + 1) throw py::value_error("TwoPortBody: len(surfaces) = len(layers) + 1 erwartet");
                 return TwoPortBody{s, core, l};
             }),
             "surfaces"_a, "core"_a, "layers"_a)
        .def_property("surfaces", [](const TwoPortBody& b) { return b.surfaces; },
                      [](TwoPortBody& b, std::vector<TriangleMesh> s) { b.surfaces = std::move(s); })
        .def_readwrite("core", &TwoPortBody::core)
        .def_property("layers", [](const TwoPortBody& b) { return b.layers; }, [](TwoPortBody& b, std::vector<Coating> l) { b.layers = std::move(l); });

    py::class_<TwoPortLayerProblem>(m, "TwoPortLayerProblem", "Schichten als Zweitor (S-Matrix-Formulierung), beliebige Dicke")
        .def(py::init([](const TriangleMesh& core_s, const TriangleMesh& outer_s, const Medium& core, const Medium& layer, real d, real omega,
                         const Medium& outer, const HMatrixParams& hp, const EntryParams& ep, const TwoPortOptions& opt) {
                 return nogil([&] { return std::make_unique<TwoPortLayerProblem>(core_s, outer_s, core, layer, d, omega, outer, hp, ep, opt); });
             }),
             "core_surface"_a, "outer_surface"_a, "core"_a, "layer"_a, "d"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = HMatrixParams{},
             "entries"_a = EntryParams{}, "options"_a = TwoPortOptions{}, "eine Schicht")
        .def(py::init([](const std::vector<TriangleMesh>& s, const Medium& core, const std::vector<Coating>& l, real omega, const Medium& outer,
                         const HMatrixParams& hp, const EntryParams& ep, const TwoPortOptions& opt) {
                 if (s.size() != l.size() + 1) throw py::value_error("len(surfaces) = len(layers) + 1 erwartet");
                 return nogil([&] { return std::make_unique<TwoPortLayerProblem>(s, core, l, omega, outer, hp, ep, opt); });
             }),
             "surfaces"_a, "core"_a, "layers"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = HMatrixParams{}, "entries"_a = EntryParams{},
             "options"_a = TwoPortOptions{}, "Mehrfachschichten (Flaechen von innen nach aussen)")
        .def(py::init([](const std::vector<TwoPortBody>& b, real omega, const Medium& outer, const HMatrixParams& hp, const EntryParams& ep,
                         const TwoPortOptions& opt) {
                 if (b.empty()) throw py::value_error("mindestens ein Koerper erwartet");
                 return nogil([&] { return std::make_unique<TwoPortLayerProblem>(b, omega, outer, hp, ep, opt); });
             }),
             "bodies"_a, "omega"_a, "outer"_a = Medium{}, "hmatrix"_a = HMatrixParams{}, "entries"_a = EntryParams{},
             "options"_a = TwoPortOptions{}, "mehrere Koerper")
        .def_static("layer_surfaces", [](const TriangleMesh& core, const std::vector<Coating>& l) {
                        return nogil([&] { return TwoPortLayerProblem::layer_surfaces(core, l); });
                    }, "core_surface"_a, "layers"_a, "Parallelflaechen Gamma_0 ... Gamma_L zu den Schichtdicken")
        .def("solve_plane_wave", &TwoPortLayerProblem::solve_plane_wave, "d"_a, "p"_a, "options"_a = SolveOptions{},
             py::call_guard<py::gil_scoped_release>())
        .def("solve_rhs", [](const TwoPortLayerProblem& P, const CArr& b, const SolveOptions& o) {
                 auto bv = to_trace(b, P.outer_mesh().size(), "b");
                 return nogil([&] { return P.solve_rhs(bv, o); });
             }, "b"_a, "options"_a = SolveOptions{}, "einfallende Spur auf outer_mesh; Aussenspur = r.h[:8 N]")
        .def("apply", [](const TwoPortLayerProblem& P, const CArr& x) {
            return apply_checked(P.size(), x, [&](const std::vector<cplx>& a, std::vector<cplx>& y) { P.apply(a, y); });
        }, "x"_a)
        .def("precondition", [](const TwoPortLayerProblem& P, const CArr& x) {
            return apply_checked(P.size(), x, [&](const std::vector<cplx>& a, std::vector<cplx>& y) { P.precondition(a, y); });
        }, "x"_a)
        .def("use_recycling", &TwoPortLayerProblem::use_recycling, "max_recycle"_a = 120)
        .def("use_gcrodr", &TwoPortLayerProblem::use_gcrodr, "k"_a = 20, "m"_a = 80)
        .def_property_readonly("recycled", &TwoPortLayerProblem::recycled)
        .def("__len__", &TwoPortLayerProblem::size)
        .def_property_readonly("unknowns", &TwoPortLayerProblem::size)
        .def_property_readonly("bodies", &TwoPortLayerProblem::bodies)
        .def("layers", &TwoPortLayerProblem::layers, "b"_a = 0, "Zahl der Schichten nach der Unterteilung")
        .def("input_layers", &TwoPortLayerProblem::input_layers, "b"_a = 0)
        .def("poles", &TwoPortLayerProblem::poles, "l"_a = 0, "b"_a = 0)
        .def_property_readonly("mean_inner_iterations", &TwoPortLayerProblem::mean_inner_iterations)
        .def_property_readonly("outer_mesh", &TwoPortLayerProblem::outer_mesh, py::return_value_policy::reference_internal,
                               "Vereinigung der Aussenflaechen");
}

}  // namespace cbem::py_bind
