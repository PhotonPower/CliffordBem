// Python-Anbindung: Anregungen, Fernfeld, Extinktion, Nahfeld, Dipole, optische Kraefte (v0.42).
#include "common.hpp"

#include "cbem/sources/chiral_incidence.hpp"
#include "cbem/sources/dipole.hpp"
#include "cbem/sources/fields.hpp"
#include "cbem/sources/incident_field.hpp"
#include "cbem/sources/near_field.hpp"
#include "cbem/sources/optical_force.hpp"

namespace cbem::py_bind {

namespace {

// Nahfeldpunkte als Woerterbuch von Arrays
py::dict near_field_to_dict(const std::vector<NearFieldPoint>& v) {
    const auto M = static_cast<py::ssize_t>(v.size());
    py::array_t<cplx> E({M, py::ssize_t(3)}), H({M, py::ssize_t(3)});
    py::array_t<bool> inside(M), too_close(M);
    py::array_t<double> enh(M), chi(M);
    cplx* e = E.mutable_data(); cplx* h = H.mutable_data();
    for (std::size_t i = 0; i < v.size(); ++i) {
        for (int c = 0; c < 3; ++c) { e[3 * i + c] = v[i].E[c]; h[3 * i + c] = v[i].H[c]; }
        inside.mutable_data()[i] = v[i].inside;
        too_close.mutable_data()[i] = v[i].too_close;
        enh.mutable_data()[i] = v[i].enhancement;
        chi.mutable_data()[i] = v[i].chirality;
    }
    return py::dict("E"_a = E, "H"_a = H, "inside"_a = inside, "too_close"_a = too_close, "enhancement"_a = enh, "chirality"_a = chi);
}

// E und H an Punkten (ein Punkt: (3,), mehrere: (M, 3))
template <class F> py::tuple eval_fields(const RArr& x, F&& f) {
    const auto pts = to_points(x);
    std::vector<CVec3> E(pts.size()), H(pts.size());
    nogil([&] { for (std::size_t i = 0; i < pts.size(); ++i) f(pts[i], E[i], H[i]); });
    if (single_point(x)) return py::make_tuple(py::cast(E[0]), py::cast(H[0]));
    return py::make_tuple(cvecs_to_numpy(E), cvecs_to_numpy(H));
}

// Feldauswerter mit eigenem Netz und eigener Spur: make_near_field_eval haelt beide nur per Referenz
struct PyNearFieldEval {
    std::shared_ptr<const TriangleMesh> mesh;
    std::shared_ptr<const std::vector<cplx>> h;
    NearFieldEval f;
    std::vector<NearFieldPoint> operator()(const std::vector<Vec3>& p) const { return f(p); }
};

py::array_t<cplx> mat33(const CVec3 (&d)[3]) {
    py::array_t<cplx> a({py::ssize_t(3), py::ssize_t(3)});
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) a.mutable_data()[3 * i + j] = d[i][j];
    return a;
}

}  // namespace

void init_sources(py::module_& m) {
    // --- ebene Wellen, Fernfeld, Extinktion -------------------------------------------------------------------------------------
    m.def("project_plane_wave", [](const TriangleMesh& me, cplx k, cplx eps, const Vec3& d, const CVec3& p, int sub) {
              return to_numpy(nogil([&] { return project_plane_wave(me, k, eps, d, p, sub); }));
          }, "mesh"_a, "k"_a, "eps"_a, "d"_a, "p"_a, "sub"_a = 2, "L2-Projektion der ebenen Welle sqrt(eps)(p + I d x p) e^(ik d.x)");
    m.def("plane_wave_trace", [](const TriangleMesh& me, const Medium& md, real omega, const Vec3& d, const CVec3& p, int sub) {
              return to_numpy(nogil([&] {
                  const PlaneWaveIncidence inc = plane_wave_incidence(md, omega, d, p);
                  return project_plane_wave(me, inc.k, md.eps, d, p, sub);
              }));
          }, "mesh"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "sub"_a = 2,
          "einfallende Spur b wie in solve_plane_wave (chirales Medium: Helizitaetswelle mit k_sigma); Streuspur = h - b");
    m.def("circular_polarization", &circular_polarization, "d"_a, "s"_a, "p = u + i s v mit (u, v, d) rechtshaendig");
    m.def("far_field", [](const TriangleMesh& me, const CArr& hs, cplx k, const Vec3& xhat, int sub) {
              auto h = to_trace(hs, me.size(), "h_scat");
              return nogil([&] { return far_field(me, h, k, xhat, sub); });
          }, "mesh"_a, "h_scat"_a, "k"_a, "xhat"_a, "sub"_a = 2, "Fernfeld-Multivektor F_inf(xhat) der Streuspur");
    m.def("far_field_E", [](const TriangleMesh& me, const CArr& hs, cplx k, cplx eps, const RArr& xhat, int sub) {
              auto h = to_trace(hs, me.size(), "h_scat");
              auto dirs = to_points(xhat);
              std::vector<CVec3> E(dirs.size());
              const cplx se = std::sqrt(eps);
              nogil([&] {
                  for (std::size_t i = 0; i < dirs.size(); ++i) {
                      const Multivector F = far_field(me, h, k, dirs[i], sub);
                      E[i] = CVec3{F.c[1] / se, F.c[2] / se, F.c[4] / se};
                  }
              });
              return cvecs_to_numpy(E);
          }, "mesh"_a, "h_scat"_a, "k"_a, "eps"_a, "directions"_a, "sub"_a = 2, "E_inf = Vektoranteil / sqrt(eps) fuer Richtungen (M, 3)");
    m.def("extinction_cross_section", [](const TriangleMesh& me, const CArr& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
              auto h = to_trace(hs, me.size(), "h_scat");
              return nogil([&] { return extinction_cross_section(me, h, k, eps, d, p); });
          }, "mesh"_a, "h_scat"_a, "k"_a, "eps"_a, "d"_a, "p"_a, "sigma_ext = 4 pi / k Im(conj(p).E_inf(d)) / |p|^2");
    m.def("forward_amplitude", [](const TriangleMesh& me, const CArr& hs, cplx k, cplx eps, const Vec3& d, const CVec3& p) {
              auto h = to_trace(hs, me.size(), "h_scat");
              return nogil([&] { return forward_amplitude(me, h, k, eps, d, p); });
          }, "mesh"_a, "h_scat"_a, "k"_a, "eps"_a, "d"_a, "p"_a, "Vorwaertsamplitude S(0) (Bohren-Huffman)");

    py::class_<PlaneWaveIncidence>(m, "PlaneWaveIncidence", "Wellenzahl und Helizitaetskanal einer ebenen Welle")
        .def_readonly("k", &PlaneWaveIncidence::k)
        .def_readonly("proj", &PlaneWaveIncidence::proj, "+1 (P+), -1 (P-), 0 (achiral)");
    m.def("plane_wave_incidence", &plane_wave_incidence, "medium"_a, "omega"_a, "d"_a, "p"_a,
          "ValueError, wenn das Medium chiral und p keine Helizitaetswelle ist");
    m.def("helicity_part", [](const CArr& h, int s) { return to_numpy(helicity_part(to_cvec(h), s)); }, "h"_a, "sigma"_a, "P_sigma h");
    m.def("extinction_in_medium", [](const TriangleMesh& me, const CArr& hs, const Medium& md, const PlaneWaveIncidence& inc, const Vec3& d,
                                     const CVec3& p) {
              auto h = to_trace(hs, me.size(), "h_scat");
              return nogil([&] { return extinction_in_medium(me, h, md, inc, d, p); });
          }, "mesh"_a, "h_scat"_a, "medium"_a, "incidence"_a, "d"_a, "p"_a);
    m.def("forward_amplitude_in_medium", [](const TriangleMesh& me, const CArr& hs, const Medium& md, const PlaneWaveIncidence& inc,
                                            const Vec3& d, const CVec3& p) {
              auto h = to_trace(hs, me.size(), "h_scat");
              return nogil([&] { return forward_amplitude_in_medium(me, h, md, inc, d, p); });
          }, "mesh"_a, "h_scat"_a, "medium"_a, "incidence"_a, "d"_a, "p"_a);

    // --- einfallende Felder -----------------------------------------------------------------------------------------------------
    py::class_<IncidentField, std::shared_ptr<IncidentField>>(m, "IncidentField", "einfallendes Feld (ebene Welle, Dipol, Strahl)")
        .def("eval", [](const IncidentField& f, const RArr& x) {
                 return eval_fields(x, [&](const Vec3& p, CVec3& E, CVec3& H) { f.eval(p, E, H); });
             }, "x"_a, "(E, H) an einem Punkt (3,) oder an Punkten (M, 3)")
        .def("project", [](const IncidentField& f, const TriangleMesh& me, const Medium& md, int sub) {
                 return to_numpy(nogil([&] { return f.project(me, md, sub); }));
             }, "mesh"_a, "medium"_a, "sub"_a = 2, "Spur sqrt(eps) E + I sqrt(mu) H auf dem Netz (rechte Seite fuer solve_rhs)")
        .def_property_readonly("reference_E2", &IncidentField::reference_E2)
        .def_property_readonly("reference_C", &IncidentField::reference_C);
    py::class_<PlaneWaveField, IncidentField, std::shared_ptr<PlaneWaveField>>(m, "PlaneWaveField")
        .def(py::init<const Medium&, real, const Vec3&, const CVec3&>(), "medium"_a, "omega"_a, "d"_a, "p"_a);
    py::class_<ZeroField, IncidentField, std::shared_ptr<ZeroField>>(m, "ZeroField", "kein einfallendes Feld (nur Streufeld)")
        .def(py::init<>());
    py::class_<DipoleField, IncidentField, std::shared_ptr<DipoleField>>(m, "DipoleField", "elektrischer (und magnetischer) Dipol")
        .def(py::init<const Medium&, real, const Vec3&, const CVec3&, const CVec3&>(), "medium"_a, "omega"_a, "r0"_a, "p"_a,
             "md"_a = CVec3{})
        .def_property_readonly("position", &DipoleField::position)
        .def_property_readonly("p", &DipoleField::p)
        .def_property_readonly("md", &DipoleField::md);
    py::class_<BeamField, IncidentField, std::shared_ptr<BeamField>>(m, "BeamField", "Strahl als Winkelspektrum, Leistung 1")
        .def_static("focused", [](const Medium& md, real omega, const Vec3& focus, real NA, real f0, const CVec3& pol, int nt, int np) {
                        return nogil([&] { return std::make_shared<BeamField>(BeamField::focused(md, omega, focus, NA, f0, pol, nt, np)); });
                    }, "medium"_a, "omega"_a, "focus"_a, "NA"_a, "f0"_a, "pupil_pol"_a, "ntheta"_a = 40, "nphi"_a = 80,
                    "fokussierter Strahl nach Richards-Wolf (NA, Fuellfaktor f0, Pupillenpolarisation)")
        .def_static("gaussian", [](const Medium& md, real omega, const Vec3& focus, real w0, const CVec3& pol, int nt, int np) {
                        return nogil([&] { return std::make_shared<BeamField>(BeamField::gaussian(md, omega, focus, w0, pol, nt, np)); });
                    }, "medium"_a, "omega"_a, "focus"_a, "w0"_a, "pol"_a, "ntheta"_a = 40, "nphi"_a = 80, "Gaussstrahl mit Taille w0")
        .def_property_readonly("power", &BeamField::power)
        .def_property_readonly("components", &BeamField::components)
        .def("poynting_flux", &BeamField::poynting_flux, "dz"_a, "half_width"_a, "ng"_a, py::call_guard<py::gil_scoped_release>(),
             "Poynting-Fluss durch die Ebene z = z_f + dz (Pruefung der Normierung)");

    // --- Nahfeld ------------------------------------------------------------------------------------------------------------------
    m.def("scattered_field", [](const TriangleMesh& me, const CArr& hs, cplx k, const RArr& pts) {
              auto h = to_trace(hs, me.size(), "h_scat");
              auto p = to_points(pts);
              return multivectors_to_numpy(nogil([&] { return scattered_field(me, h, k, p); }));
          }, "mesh"_a, "h_scat"_a, "k"_a, "points"_a, "Streufeld-Multivektoren (M, 8), ein Helizitaetskanal, direkte Summation");
    py::class_<NearFieldOperator>(m, "NearFieldOperator", "Nahfeldoperator als rechteckige H-Matrix (Punkte x Dreiecke)")
        .def(py::init([](const TriangleMesh& me, cplx k, const RArr& pts, const HMatrixParams& hp) {
                 auto p = to_points(pts);
                 return nogil([&] { return std::make_unique<NearFieldOperator>(me, k, p, hp); });
             }),
             "mesh"_a, "k"_a, "points"_a, "hmatrix"_a = HMatrixParams{}, py::keep_alive<1, 2>())
        .def("apply", [](const NearFieldOperator& N, const CArr& hs) {
                 auto h = to_cvec(hs);
                 return multivectors_to_numpy(nogil([&] { return N.apply(h); }));
             }, "h_scat"_a, "Streufeld-Multivektoren (M, 8)")
        .def_property_readonly("stats", &NearFieldOperator::stats);
    m.def("exterior_near_field", [](const TriangleMesh& outer, const CArr& h, const Medium& md, real omega, const Vec3& d, const CVec3& p,
                                    const RArr& pts, const NearFieldOptions& o) {
              auto hv = to_trace(h, outer.size(), "h");
              auto x = to_points(pts);
              return near_field_to_dict(nogil([&] { return exterior_near_field(outer, hv, md, omega, d, p, x, o); }));
          }, "outer"_a, "h"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "points"_a, "options"_a = NearFieldOptions{}, R"doc(
Gesamtes Nahfeld bei Anregung durch eine ebene Welle; h = Gesamtspur auf der Flaeche zum Aussenraum.
Gibt ein dict mit E, H (M, 3), inside, too_close (M,), enhancement = |E|^2/|E0|^2 und chirality (M,) zurueck.
)doc");
    m.def("exterior_near_field", [](const TriangleMesh& outer, const CArr& h, const CArr& b, const Medium& md, real omega,
                                    const IncidentField& inc, const RArr& pts, const NearFieldOptions& o) {
              auto hv = to_trace(h, outer.size(), "h");
              auto bv = to_trace(b, outer.size(), "b");
              auto x = to_points(pts);
              return near_field_to_dict(nogil([&] { return exterior_near_field(outer, hv, bv, md, omega, inc, x, o); }));
          }, "outer"_a, "h"_a, "b"_a, "medium"_a, "omega"_a, "incident"_a, "points"_a, "options"_a = NearFieldOptions{},
          "allgemeine Anregung: b = Projektion des einfallenden Feldes (incident.project oder project_dipole)");

    py::class_<PyNearFieldEval>(m, "NearFieldEvaluator", R"doc(
Feldauswerter fuer Kraefte und Gradienten: haelt eigene Kopien von Netz und Spur. Aufruf mit Punkten (M, 3)
liefert dasselbe dict wie exterior_near_field.
)doc")
        .def("__call__", [](const PyNearFieldEval& e, const RArr& pts) {
            auto x = to_points(pts);
            return near_field_to_dict(nogil([&] { return e(x); }));
        }, "points"_a);
    m.def("near_field_evaluator", [](const TriangleMesh& outer, const CArr& h, const CArr& b, const Medium& md, real omega,
                                     std::shared_ptr<IncidentField> inc, const NearFieldOptions& o) {
              PyNearFieldEval e;
              e.mesh = std::make_shared<const TriangleMesh>(outer);
              e.h = std::make_shared<const std::vector<cplx>>(to_trace(h, outer.size(), "h"));
              e.f = make_near_field_eval(*e.mesh, *e.h, to_trace(b, outer.size(), "b"), md, omega, inc, o);
              return e;
          }, "outer"_a, "h"_a, "b"_a, "medium"_a, "omega"_a, "incident"_a, "options"_a = NearFieldOptions{},
          "Feldauswerter fuer beliebige Anregung (Strahl, Dipol, ebene Welle)");
    m.def("plane_wave_evaluator", [](const TriangleMesh& outer, const CArr& h, const Medium& md, real omega, const Vec3& d, const CVec3& p,
                                     const NearFieldOptions& o) {
              PyNearFieldEval e;
              e.mesh = std::make_shared<const TriangleMesh>(outer);
              e.h = std::make_shared<const std::vector<cplx>>(to_trace(h, outer.size(), "h"));
              e.f = make_plane_wave_eval(*e.mesh, *e.h, md, omega, d, p, o);
              return e;
          }, "outer"_a, "h"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "options"_a = NearFieldOptions{}, "Feldauswerter fuer eine ebene Welle");

    // --- Dipole ---------------------------------------------------------------------------------------------------------------------
    m.def("dipole_field", [](const Medium& md, real omega, const Vec3& r0, const CVec3& p, const RArr& x, const CVec3& mdip) {
              return eval_fields(x, [&](const Vec3& q, CVec3& E, CVec3& H) { dipole_field(md, omega, r0, p, q, E, H, mdip); });
          }, "medium"_a, "omega"_a, "r0"_a, "p"_a, "x"_a, "md"_a = CVec3{}, "(E, H) des Dipols an einem Punkt oder an Punkten (M, 3)");
    m.def("project_dipole", [](const TriangleMesh& me, const Medium& md, real omega, const Vec3& r0, const CVec3& p, int max_sub,
                               const CVec3& mdip) {
              return to_numpy(nogil([&] { return project_dipole(me, md, omega, r0, p, max_sub, mdip); }));
          }, "mesh"_a, "medium"_a, "omega"_a, "r0"_a, "p"_a, "max_sub"_a = 24, "md"_a = CVec3{},
          "Spur des Dipolfeldes, nahe Dreiecke adaptiv unterteilt");
    py::class_<DipoleRates>(m, "DipoleRates", "Zerfallsraten relativ zum freien Dipol")
        .def(py::init<>())
        .def_readwrite("total", &DipoleRates::total)
        .def_readwrite("radiative", &DipoleRates::radiative)
        .def_readwrite("nonradiative", &DipoleRates::nonradiative)
        .def_readwrite("distance", &DipoleRates::distance)
        .def_readwrite("too_close", &DipoleRates::too_close)
        .def_readwrite("rad_plus", &DipoleRates::rad_plus)
        .def_readwrite("rad_minus", &DipoleRates::rad_minus)
        .def_readwrite("glum", &DipoleRates::glum)
        .def_readwrite("glum_free", &DipoleRates::glum_free)
        .def_property_readonly("quantum_yield", [](const DipoleRates& r) { return r.total > 0 ? r.radiative / r.total : 0.0; },
                               "radiative / total (intrinsische Quantenausbeute 1)")
        .def("__repr__", [](const DipoleRates& r) {
            auto f = [](real v) { return py::repr(py::float_(v)).cast<std::string>(); };
            return "DipoleRates(total=" + f(r.total) + ", radiative=" + f(r.radiative) + ", nonradiative=" + f(r.nonradiative) +
                   ", distance=" + f(r.distance) + ", glum=" + f(r.glum) + ")";
        });
    m.def("dipole_rates", [](const TriangleMesh& outer, const CArr& h, const CArr& b, const Medium& md, real omega, const Vec3& r0,
                             const CVec3& p, int ntheta, const CVec3& mdip) {
              auto hv = to_trace(h, outer.size(), "h");
              auto bv = to_trace(b, outer.size(), "b");
              return nogil([&] { return dipole_rates(outer, hv, bv, md, omega, r0, p, ntheta, mdip); });
          }, "outer"_a, "h"_a, "b"_a, "medium"_a, "omega"_a, "r0"_a, "p"_a, "ntheta"_a = 40, "md"_a = CVec3{});
    m.def("helicity_projector_sign", &helicity_projector_sign, "s"_a);
    m.def("fluorescence_enhancement", [](const std::vector<real>& exc, const std::vector<DipoleRates>& rates, real q0) {
              if (exc.size() != 3 || rates.size() != 3) throw py::value_error("je drei Werte (x, y, z) erwartet");
              const real e[3] = {exc[0], exc[1], exc[2]};
              const DipoleRates r[3] = {rates[0], rates[1], rates[2]};
              return fluorescence_enhancement(e, r, q0);
          }, "exc"_a, "rates"_a, "q0"_a, "F/F0 eines fest, aber zufaellig orientierten Emitters");

    // --- optische Kraefte ---------------------------------------------------------------------------------------------------------
    m.def("stress_dot_normal", &stress_dot_normal, "E"_a, "H"_a, "n"_a, "medium"_a, "<T>.n (Minkowski, auch chiral)");
    m.def("force_from_traces", [](const TriangleMesh& outer, const CArr& h, const Medium& md, const std::vector<std::size_t>& bb) {
              auto hv = to_trace(h, outer.size(), "h");
              return points_to_numpy(nogil([&] { return force_from_traces(outer, hv, md, bb); }));
          }, "outer"_a, "h"_a, "medium"_a, "body_begin"_a, "Kraft je Koerper (B, 3) aus den Randspuren");
    m.def("force_from_far_field", [](const TriangleMesh& outer, const CArr& h, const Medium& md, real omega, const Vec3& d, const CVec3& p,
                                     real sext, int nt) {
              auto hv = to_trace(h, outer.size(), "h");
              return nogil([&] { return force_from_far_field(outer, hv, md, omega, d, p, sext, nt); });
          }, "outer"_a, "h"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "sigma_ext"_a, "ntheta"_a = 32, "Kraft aus der Impulsbilanz im Fernfeld");
    m.def("force_on_sphere", [](const PyNearFieldEval& e, const Medium& md, const Vec3& c, real R, int nt) {
              return nogil([&] { return force_on_sphere(e.f, md, c, R, nt); });
          }, "evaluator"_a, "medium"_a, "center"_a, "R"_a, "ntheta"_a = 32, "Kraft auf alles innerhalb der Kugel (Spannungstensor)");
    m.def("force_on_sphere", [](const TriangleMesh& outer, const CArr& h, const Medium& md, real omega, const Vec3& d, const CVec3& p,
                                const Vec3& c, real R, int nt, const NearFieldOptions& o) {
              auto hv = to_trace(h, outer.size(), "h");
              return nogil([&] { return force_on_sphere(outer, hv, md, omega, d, p, c, R, nt, o); });
          }, "outer"_a, "h"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "center"_a, "R"_a, "ntheta"_a = 32, "options"_a = NearFieldOptions{});
    m.def("force_on_offset", [](const PyNearFieldEval& e, const Medium& md, const TriangleMesh& body, real delta) {
              return nogil([&] { return force_on_offset(e.f, md, body, delta); });
          }, "evaluator"_a, "medium"_a, "body"_a, "delta"_a, "Kraft auf den Koerper ueber die Parallelflaeche im Abstand delta");
    m.def("force_on_offset", [](const TriangleMesh& outer, const CArr& h, const Medium& md, real omega, const Vec3& d, const CVec3& p,
                                const TriangleMesh& body, real delta, const NearFieldOptions& o) {
              auto hv = to_trace(h, outer.size(), "h");
              return nogil([&] { return force_on_offset(outer, hv, md, omega, d, p, body, delta, o); });
          }, "outer"_a, "h"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "body"_a, "delta"_a, "options"_a = NearFieldOptions{});
    m.def("emitter_force", [](const TriangleMesh& outer, const CArr& h, const CArr& b, const Medium& md, real omega, const DipoleField& src,
                              real delta) {
              auto hv = to_trace(h, outer.size(), "h");
              auto bv = to_trace(b, outer.size(), "b");
              return nogil([&] { return emitter_force(outer, hv, bv, md, omega, src, delta); });
          }, "outer"_a, "h"_a, "b"_a, "medium"_a, "omega"_a, "source"_a, "delta"_a, "Kraft auf den Emitter (nur Streufeld)");
    m.def("radiated_momentum", [](const TriangleMesh& outer, const CArr& h, const CArr& b, const Medium& md, real omega, const DipoleField& src,
                                  int nt) {
              auto hv = to_trace(h, outer.size(), "h");
              auto bv = to_trace(b, outer.size(), "b");
              return nogil([&] { return radiated_momentum(outer, hv, bv, md, omega, src, nt); });
          }, "outer"_a, "h"_a, "b"_a, "medium"_a, "omega"_a, "source"_a, "ntheta"_a = 32, "abgestrahlter Impuls (Impulsbilanz)");

    py::class_<DipolePolarizability>(m, "DipolePolarizability", "Polarisierbarkeiten A_e, A_m, A_c (duale Groessen)")
        .def(py::init([](cplx Ae, cplx Am, cplx Ac) { return DipolePolarizability{Ae, Am, Ac}; }), "Ae"_a = cplx(0), "Am"_a = cplx(0),
             "Ac"_a = cplx(0))
        .def_readwrite("Ae", &DipolePolarizability::Ae)
        .def_readwrite("Am", &DipolePolarizability::Am)
        .def_readwrite("Ac", &DipolePolarizability::Ac);
    m.def("polarizability_from_mie", &polarizability_from_mie, "a1"_a, "b1"_a, "k"_a, "A_e = 6 pi i a_1 / k^3, A_m = 6 pi i b_1 / k^3");
    py::class_<FieldGradient>(m, "FieldGradient", "E, H und Gradienten dE[i, j] = d_i E_j")
        .def(py::init<>())
        .def_readonly("E", &FieldGradient::E)
        .def_readonly("H", &FieldGradient::H)
        .def_property_readonly("dE", [](const FieldGradient& g) { return mat33(g.dE); })
        .def_property_readonly("dH", [](const FieldGradient& g) { return mat33(g.dH); });
    m.def("fields_with_gradients", [](const PyNearFieldEval& e, const RArr& x, real delta) {
              auto p = to_points(x);
              return nogil([&] { return fields_with_gradients(e.f, p, delta); });
          }, "evaluator"_a, "points"_a, "delta"_a, "Felder und Gradienten (zentrale Differenzen)");
    m.def("fields_with_gradients", [](const TriangleMesh& outer, const CArr& h, const Medium& md, real omega, const Vec3& d, const CVec3& p,
                                      const RArr& x, real delta, const NearFieldOptions& o) {
              auto hv = to_trace(h, outer.size(), "h");
              auto q = to_points(x);
              return nogil([&] { return fields_with_gradients(outer, hv, md, omega, d, p, q, delta, o); });
          }, "outer"_a, "h"_a, "medium"_a, "omega"_a, "d"_a, "p"_a, "points"_a, "delta"_a, "options"_a = NearFieldOptions{});
    m.def("dipole_particle_force", &dipole_particle_force, "g"_a, "alpha"_a, "medium"_a, "omega"_a,
          "Kraft auf ein kleines (auch chirales) Teilchen in Dipolnaeherung");
}

}  // namespace cbem::py_bind
