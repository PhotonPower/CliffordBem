// Python-Anbindung: Grundtypen, Clifford-Algebra, Medien, Materialien, Parameterstrukturen (v0.42).
#include "common.hpp"

#include <sstream>

#include "cbem/assembly/kernel_entries.hpp"
#include "cbem/core/materials.hpp"
#include "cbem/hmatrix/hmatrix.hpp"
#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/problems/thin_layer_problem.hpp"
#include "cbem/problems/twoport_layer_problem.hpp"
#include "cbem/solvers/hodlr.hpp"
#include "cbem/sources/near_field.hpp"

namespace cbem::py_bind {

namespace {

const char* const kBladeNames[8] = {"1", "e1", "e2", "e12", "e3", "e13", "e23", "e123"};

int blade_index(const std::string& s) {
    for (int b = 0; b < 8; ++b) if (s == kBladeNames[b]) return b;
    if (s == "I") return 7;
    throw py::key_error("unbekanntes Blade '" + s + "' (erlaubt: 1, e1, e2, e12, e3, e13, e23, e123, I)");
}

std::string complex_repr(cplx z) { return py::repr(py::cast(z)).cast<std::string>(); }

std::string mv_repr(const Multivector& M) {
    std::string s;
    for (int b = 0; b < 8; ++b) {
        if (M.c[b] == cplx(0)) continue;
        if (!s.empty()) s += " + ";
        const cplx z = M.c[b];
        std::string v = z.imag() == 0 ? py::repr(py::float_(z.real())).cast<std::string>() : complex_repr(z);
        s += b == 0 ? v : v + "*" + kBladeNames[b];
    }
    return "Multivector(" + (s.empty() ? std::string("0") : s) + ")";
}

Multivector grade_part(const Multivector& M, int g) {
    Multivector r;
    for (int b = 0; b < 8; ++b) if (blade_grade(b) == g) r.c[b] = M.c[b];
    return r;
}

std::string medium_repr(const Medium& m) {
    return "Medium(eps=" + complex_repr(m.eps) + ", mu=" + complex_repr(m.mu) + ", chi=" + complex_repr(m.chi) + ")";
}

}  // namespace

void init_core(py::module_& m) {
    // --- Multivektoren ------------------------------------------------------------------------------------------------------
    py::class_<Multivector>(m, "Multivector", R"doc(
Element der komplexifizierten Clifford-Algebra Cl3(C).

Koeffizienten in der Blade-Basis (Bitmasken) 1, e1, e2, e12, e3, e13, e23, e123. Die imaginaere Einheit i
(komplexe Koeffizienten) ist vom Pseudoskalar I = e123 verschieden. Multivektoren sind unveraenderlich;
Operationen liefern neue Objekte. ``*`` ist das geometrische Produkt.

Beispiel::

    e1, e2 = Multivector.blade(1), Multivector.blade(2)
    (e1 * e2)["e12"]          # 1
    I = Multivector.blade(7); (I * I)["1"]   # -1
)doc")
        .def(py::init<>(), "Null-Multivektor")
        .def(py::init([](const CArr& c) {
                 if (c.size() != 8) throw py::value_error("Multivector: genau 8 Koeffizienten erwartet");
                 Multivector M; std::copy(c.data(), c.data() + 8, M.c.begin()); return M;
             }),
             "coefficients"_a, "aus 8 Koeffizienten (Blade-Reihenfolge 1, e1, e2, e12, e3, e13, e23, e123)")
        .def_static("blade", [](py::object b, cplx v) {
                        const int idx = py::isinstance<py::str>(b) ? blade_index(b.cast<std::string>()) : b.cast<int>();
                        if (idx < 0 || idx > 7) throw py::index_error("Blade-Index 0..7");
                        return Multivector::blade(idx, v);
                    },
                    "blade"_a, "value"_a = cplx(1.0), "Blade als Index (Bitmaske 0..7) oder Name ('e1', 'e12', 'I', ...)")
        .def_static("scalar", [](cplx v) { return Multivector::blade(0, v); }, "value"_a)
        .def_static("vector", [](const CVec3& v) { return Multivector::vector(v); }, "v"_a, "Vektor v1 e1 + v2 e2 + v3 e3")
        .def_property_readonly("c", [](const Multivector& M) {
                                   py::array_t<cplx> a(8); std::copy(M.c.begin(), M.c.end(), a.mutable_data()); return a; },
                               "Koeffizienten (Kopie)")
        .def("__getitem__", [](const Multivector& M, py::object k) {
                 const int idx = py::isinstance<py::str>(k) ? blade_index(k.cast<std::string>()) : k.cast<int>();
                 if (idx < 0 || idx > 7) throw py::index_error("Blade-Index 0..7");
                 return M.c[idx];
             })
        .def("__len__", [](const Multivector&) { return 8; })
        .def("__add__", [](const Multivector& a, const Multivector& b) { return a + b; })
        .def("__add__", [](const Multivector& a, cplx s) { return a + Multivector::blade(0, s); })
        .def("__radd__", [](const Multivector& a, cplx s) { return a + Multivector::blade(0, s); })
        .def("__sub__", [](const Multivector& a, const Multivector& b) { return a - b; })
        .def("__sub__", [](const Multivector& a, cplx s) { return a - Multivector::blade(0, s); })
        .def("__rsub__", [](const Multivector& a, cplx s) { return Multivector::blade(0, s) - a; })
        .def("__neg__", [](const Multivector& a) { return a * cplx(-1.0); })
        .def("__mul__", [](const Multivector& a, const Multivector& b) { return a * b; })
        .def("__mul__", [](const Multivector& a, cplx s) { return a * s; })
        .def("__rmul__", [](const Multivector& a, cplx s) { return a * s; })
        .def("__truediv__", [](const Multivector& a, cplx s) { return a * (cplx(1.0) / s); })
        .def("__eq__", [](const Multivector& a, const Multivector& b) { return a.c == b.c; })
        .def("__repr__", &mv_repr)
        .def("reverse", &Multivector::reverse, "Reversion (Grad 2 und 3 wechseln das Vorzeichen)")
        .def("involute", &Multivector::involute, "Grad-Involution (ungerade Grade wechseln das Vorzeichen)")
        .def("grade", &grade_part, "k"_a, "Anteil vom Grad k")
        .def("scalar_part", [](const Multivector& M) { return M.c[0]; })
        .def("vector_part", [](const Multivector& M) { return CVec3{M.c[1], M.c[2], M.c[4]}; }, "Vektoranteil (x, y, z)")
        .def("pseudoscalar_part", [](const Multivector& M) { return M.c[7]; })
        .def("left_matrix", [](const Multivector& M) { return mat8_to_numpy(M.left_matrix()); },
             "Linksmultiplikation L(M) als (8, 8)-Matrix: (M X).c = L(M) @ X.c")
        .def("inverse", [](const Multivector& M) {
                 bool ok = false; Multivector r = mv_inverse(M, &ok);
                 if (!ok) throw py::value_error("Multivektor ist ein Nullteiler (nicht invertierbar)");
                 return r;
             })
        .def("norm2", [](const Multivector& M) { return mv_norm2(M); }, "Summe der Betragsquadrate der Koeffizienten")
        .def("__copy__", [](const Multivector& M) { return M; })
        .def("__deepcopy__", [](const Multivector& M, py::dict) { return M; }, "memo"_a)
        .def(py::pickle([](const Multivector& M) { py::array_t<cplx> a(8); std::copy(M.c.begin(), M.c.end(), a.mutable_data()); return py::make_tuple(a); },
                        [](py::tuple t) { CArr a = t[0].cast<CArr>(); Multivector M; std::copy(a.data(), a.data() + 8, M.c.begin()); return M; }));
    m.def("blade_grade", &blade_grade, "b"_a, "Grad des Blades mit Bitmaske b");
    m.def("blade_sign", &blade_sign, "a"_a, "b"_a, "Vorzeichen von e_a e_b = s e_(a xor b)");
    m.attr("BLADE_NAMES") = py::make_tuple("1", "e1", "e2", "e12", "e3", "e13", "e23", "e123");

    // --- Medien und Materialien ---------------------------------------------------------------------------------------------
    py::class_<Medium>(m, "Medium", R"doc(
Pasteur-Medium: D = eps E + i chi H, B = mu H - i chi E (chi = 0: achiral). Einheiten mit eps0 = mu0 = 1,
Zeitkonvention exp(-i omega t), Absorption Im eps > 0.
)doc")
        .def(py::init([](cplx eps, cplx mu, cplx chi) { return Medium{eps, mu, chi}; }), "eps"_a = cplx(1.0), "mu"_a = cplx(1.0),
             "chi"_a = cplx(0.0))
        .def_readwrite("eps", &Medium::eps)
        .def_readwrite("mu", &Medium::mu)
        .def_readwrite("chi", &Medium::chi)
        .def("k", &Medium::k, "omega"_a, "helicity"_a = 0,
             "Wellenzahl; helicity = +-1: k_pm = omega (sqrt(eps mu) -+ chi), 0: achiral")
        .def_property_readonly("chiral", [](const Medium& me) { return me.chi != cplx(0); })
        .def("__repr__", &medium_repr)
        .def("__copy__", [](const Medium& me) { return me; })
        .def("__deepcopy__", [](const Medium& me, py::dict) { return me; }, "memo"_a)
        .def(py::pickle([](const Medium& me) { return py::make_tuple(me.eps, me.mu, me.chi); },
                        [](py::tuple t) { return Medium{t[0].cast<cplx>(), t[1].cast<cplx>(), t[2].cast<cplx>()}; }));

    py::class_<Material, std::shared_ptr<Material>>(m, "Material", "Materialmodell eps(lambda_nm)")
        .def("eps", &Material::eps, "lambda_nm"_a, "Permittivitaet bei der Vakuumwellenlaenge lambda (nm)")
        .def("eps", [](const Material& mat, const RArr& lam) {
                 py::array_t<cplx> r(lam.size());
                 for (py::ssize_t i = 0; i < lam.size(); ++i) r.mutable_data()[i] = mat.eps(lam.data()[i]);
                 return r;
             }, "lambda_nm"_a, "vektorisiert")
        .def("medium", [](const Material& mat, real lam, cplx chi) { return Medium{mat.eps(lam), 1.0, chi}; }, "lambda_nm"_a,
             "chi"_a = cplx(0.0), "Medium bei der Wellenlaenge lambda (nicht magnetisch)")
        .def_property_readonly("name", &Material::name)
        .def("__repr__", [](const Material& mat) { return "Material('" + mat.name() + "')"; });
    py::class_<ConstantMaterial, Material, std::shared_ptr<ConstantMaterial>>(m, "ConstantMaterial")
        .def(py::init<cplx>(), "eps"_a);
    py::class_<TabulatedMaterial, Material, std::shared_ptr<TabulatedMaterial>>(m, "TabulatedMaterial",
                                                                               "tabellierte (n, k) ueber lambda in nm, lineare Interpolation")
        .def(py::init([](std::string name, const RArr& l, const RArr& n, const RArr& k) {
                 if (l.size() != n.size() || l.size() != k.size() || l.size() < 2)
                     throw py::value_error("TabulatedMaterial: lambda, n, k gleicher Laenge (>= 2)");
                 return std::make_shared<TabulatedMaterial>(std::move(name), to_rvec(l), to_rvec(n), to_rvec(k));
             }),
             "name"_a, "lambda_nm"_a, "n"_a, "k"_a)
        .def_static("from_yaml", &TabulatedMaterial::from_yaml, "path"_a, "refractiveindex.info-YAML ('tabulated nk', lambda in um)")
        .def_property_readonly("lambda_min", &TabulatedMaterial::lambda_min)
        .def_property_readonly("lambda_max", &TabulatedMaterial::lambda_max);
    m.def("_make_material", &make_material, "spec"_a, "data_dir"_a);

    // --- Parameterstrukturen ------------------------------------------------------------------------------------------------
    py::enum_<AcaMode>(m, "AcaMode")
        .value("Joint", AcaMode::Joint).value("Componentwise", AcaMode::Componentwise).value("Multivector", AcaMode::Multivector);
    py::enum_<ThinLayerModel>(m, "ThinLayerModel")
        .value("Jump1", ThinLayerModel::Jump1).value("Dirac1", ThinLayerModel::Dirac1)
        .value("Dirac2", ThinLayerModel::Dirac2).value("Dirac2Fit", ThinLayerModel::Dirac2Fit);
    py::enum_<Adjacency>(m, "Adjacency")
        .value("None_", Adjacency::None).value("Vertex", Adjacency::Vertex).value("Edge", Adjacency::Edge)
        .value("Coincident", Adjacency::Coincident);

    ParamClass<HMatrixParams>(m, "HMatrixParams", "Parameter der H-Matrix-Kompression")
        .field("eps", &HMatrixParams::eps, "relative ACA-Toleranz je Block")
        .field("eta", &HMatrixParams::eta, "Zulaessigkeit min(diam) <= eta dist")
        .field("leaf", &HMatrixParams::leaf, "maximale Blattgroesse")
        .field("sep_factor", &HMatrixParams::sep_factor, "zusaetzlich dist > sep_factor * h_max; 0 = aus")
        .field("exact_in_lowrank", &HMatrixParams::exact_in_lowrank, "ACA mit exakten Eintraegen")
        .field("max_kdiam", &HMatrixParams::max_kdiam, "|k| diam <= max_kdiam")
        .field("mode", &HMatrixParams::mode, "AcaMode")
        .field("aca_plus", &HMatrixParams::aca_plus, "ACA+ statt teilpivotisierter ACA")
        .field("single_precision", &HMatrixParams::single_precision, "Eintraege als complex<float> speichern (halber Speicher)");
    ParamClass<EntryParams>(m, "EntryParams", "Parameter der Eintragsauswertung (Nahfeldquadratur)")
        .field("near_factor", &EntryParams::near_factor, "Nahpaar, wenn Schwerpunktabstand < near_factor * max(h_i, h_j)")
        .field("near_subdivision", &EntryParams::near_subdivision, "feste aeussere Regel: sub^2 * 7 Punkte")
        .field("sauter_schwab", &EntryParams::sauter_schwab, "benachbarte Paare mit Sauter-Schwab")
        .field("ss_order", &EntryParams::ss_order, "Gauss-Punkte je Richtung (Sauter-Schwab)")
        .field("cache_near", &EntryParams::cache_near, "Nahpaare vorab berechnen")
        .field("adaptive_outer", &EntryParams::adaptive_outer, "adaptive aeussere Regel")
        .field("adapt_ratio", &EntryParams::adapt_ratio, "Umkreisradius < adapt_ratio * Abstand")
        .field("adapt_depth", &EntryParams::adapt_depth, "maximale Halbierungstiefe")
        .field("adapt_to_boundary", &EntryParams::adapt_to_boundary, "Randabstand fuer parallele Flaechen")
        .field("sa_order", &EntryParams::sa_order, "halbanalytische Regel: Gauss-Punkte je Richtung")
        .field("ss_max_aspect", &EntryParams::ss_max_aspect, "Sauter-Schwab nur bis zu diesem Seitenverhaeltnis");
    m.def("layered_entry_params", &layered_entry_params, "Voreinstellung fuer geschichtete Koerper (adapt_to_boundary)");
    ParamClass<SolveOptions>(m, "SolveOptions", "GMRES-Optionen")
        .field("tol", &SolveOptions::tol, "relative Residuumstoleranz")
        .field("restart", &SolveOptions::restart, "Neustartlaenge")
        .field("max_iter", &SolveOptions::max_iter, "hoechstens so viele Iterationen");
    ParamClass<HodlrParams>(m, "HodlrParams", "Parameter der HODLR-Faktorisierung")
        .field("leaf", &HodlrParams::leaf, "Blattgroesse (Dreiecke)")
        .field("eps", &HodlrParams::eps, "Toleranz der Nebendiagonalbloecke");
    ParamClass<NearFieldOptions>(m, "NearFieldOptions", "Optionen der Nahfeldauswertung")
        .field("hmatrix_min_points", &NearFieldOptions::hmatrix_min_points, "H-Matrix ab so vielen Punkten (0: immer)")
        .field("eps", &NearFieldOptions::eps, "ACA-Toleranz");
    ParamClass<TwoPortOptions>(m, "TwoPortOptions", "Optionen der Zweitor-Formulierung")
        .field("poles", &TwoPortOptions::poles, "Pole je Schicht; -1: automatisch")
        .field("inner_tol", &TwoPortOptions::inner_tol, "GMRES-Toleranz der Resolventen")
        .field("curvature", &TwoPortOptions::curvature, "B in der Schichtmitte (Kruemmung)")
        .field("split", &TwoPortOptions::split, "Unterteilung dicker Schichten (0: aus)");

    py::class_<HStats>(m, "HStats", "Statistik einer H-Matrix")
        .def_readonly("n_dense", &HStats::n_dense)
        .def_readonly("n_lowrank", &HStats::n_lowrank)
        .def_readonly("entries_dense", &HStats::entries_dense)
        .def_readonly("entries_lowrank", &HStats::entries_lowrank)
        .def_readonly("mean_rank", &HStats::mean_rank)
        .def_readonly("max_rank", &HStats::max_rank)
        .def_readonly("seconds", &HStats::seconds)
        .def_property_readonly("bytes", &HStats::bytes)
        .def_readonly("entry_bytes", &HStats::entry_bytes)
        .def("__repr__", [](const HStats& s) {
            std::ostringstream o;
            o << "HStats(dense=" << s.n_dense << ", lowrank=" << s.n_lowrank << ", mean_rank=" << s.mean_rank << ", max_rank=" << s.max_rank
              << ", MB=" << s.bytes() / 1048576.0 << ", seconds=" << s.seconds << ")";
            return o.str();
        });

    // --- Transmissionsabbildung und Projektoren -------------------------------------------------------------------------------
    m.def("transmission_map", [](const Vec3& n, const Medium& inner, const Medium& outer) { return mat8_to_numpy(transmission_map(n, inner, outer)); },
          "n"_a, "inner"_a, "outer"_a, "J als (8, 8)-Matrix zur Normalen n (Medium 1 innen, Medium 2 aussen)");
    m.def("helicity_projector", [](int s) { return mat8_to_numpy(helicity_projector(s)); }, "sign"_a,
          "zentraler Helizitaetsprojektor P_pm = (1 +- iI)/2 als (8, 8)-Matrix");
}

}  // namespace cbem::py_bind
