// Python-Anbindung: Dreiecksnetze, Netzgeneratoren, Gmsh, Parallelflaechen, Richtungsquadraturen, Gruppierung (v0.42).
#include "common.hpp"

#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/geometry/mesh.hpp"
#include "cbem/solvers/block_preconditioner.hpp"
#include "cbem/sources/orientation.hpp"

namespace cbem::py_bind {

namespace {

TriangleMesh mesh_from_arrays(const RArr& P, const py::array_t<long long, py::array::c_style | py::array::forcecast>& T, bool orient,
                              const Vec3& center) {
    if (P.ndim() != 2 || P.shape(1) != 3) throw py::value_error("points: Form (n, 3) erwartet");
    if (T.ndim() != 2 || T.shape(1) != 3) throw py::value_error("triangles: Form (N, 3) erwartet");
    TriangleMesh m;
    m.P = to_points(P);
    const auto n = static_cast<long long>(m.P.size());
    m.T.resize(static_cast<std::size_t>(T.shape(0)));
    const long long* t = T.data();
    for (std::size_t i = 0; i < m.T.size(); ++i)
        for (int c = 0; c < 3; ++c) {
            const long long v = t[3 * i + c];
            if (v < 0 || v >= n) throw py::index_error("triangles: Knotenindex " + std::to_string(v) + " ausserhalb 0.." + std::to_string(n - 1));
            m.T[i][c] = static_cast<int>(v);
        }
    if (orient) m.orient_outward(center);
    m.compute_geometry();
    return m;
}

py::array_t<int> triangles_to_numpy(const TriangleMesh& m) {
    py::array_t<int> a({static_cast<py::ssize_t>(m.T.size()), py::ssize_t(3)});
    int* d = a.mutable_data();
    for (std::size_t i = 0; i < m.T.size(); ++i) for (int c = 0; c < 3; ++c) d[3 * i + c] = m.T[i][c];
    return a;
}

}  // namespace

void init_geometry(py::module_& m) {
    py::class_<TriangleMesh>(m, "TriangleMesh", R"doc(
Dreiecksnetz aus ebenen Dreiecken (gegen den Uhrzeigersinn um die Aussennormale).

Die Arrays (points, triangles, normals, centroids, areas, hmax) sind Kopien; ein Netz wird nicht an Ort und
Stelle veraendert, sondern neu erzeugt (z. B. ``TriangleMesh(2 * m.points, m.triangles)`` oder ``translated``).
)doc")
        .def(py::init(&mesh_from_arrays), "points"_a, "triangles"_a, "orient_outward"_a = false, "center"_a = Vec3(0, 0, 0),
             "Netz aus Knoten (n, 3) und Dreiecken (N, 3); orient_outward: fuer sternfoermige Gebiete um center nach aussen drehen")
        .def_property_readonly("points", [](const TriangleMesh& me) { return points_to_numpy(me.P); }, "Knoten (n, 3)")
        .def_property_readonly("triangles", &triangles_to_numpy, "Dreiecke (N, 3), Knotenindizes")
        .def_property_readonly("normals", [](const TriangleMesh& me) { return points_to_numpy(me.normal); }, "Aussennormalen (N, 3)")
        .def_property_readonly("centroids", [](const TriangleMesh& me) { return points_to_numpy(me.centroid); }, "Schwerpunkte (N, 3)")
        .def_property_readonly("areas", [](const TriangleMesh& me) { return to_numpy(me.area); }, "Flaecheninhalte (N,)")
        .def_property_readonly("hmax", [](const TriangleMesh& me) { return to_numpy(me.hmax); }, "laengste Kante je Dreieck (N,)")
        .def_property_readonly("n_points", [](const TriangleMesh& me) { return me.P.size(); })
        .def("__len__", &TriangleMesh::size, "Zahl der Dreiecke")
        .def_property_readonly("unknowns", [](const TriangleMesh& me) { return 8 * me.size(); }, "8 N (Laenge einer Spur)")
        .def("vertices", [](const TriangleMesh& me, std::size_t t) {
                 if (t >= me.size()) throw py::index_error("Dreiecksindex");
                 auto v = me.vertices(t); return points_to_numpy({v[0], v[1], v[2]});
             }, "t"_a, "Ecken des Dreiecks t als (3, 3)")
        .def("copy", [](const TriangleMesh& me) { return me; })
        .def("__copy__", [](const TriangleMesh& me) { return me; })
        .def("__deepcopy__", [](const TriangleMesh& me, py::dict) { return me; }, "memo"_a)
        .def("__repr__", [](const TriangleMesh& me) {
            return "TriangleMesh(" + std::to_string(me.P.size()) + " Knoten, " + std::to_string(me.size()) + " Dreiecke)";
        })
        .def(py::pickle([](const TriangleMesh& me) { return py::make_tuple(points_to_numpy(me.P), triangles_to_numpy(me)); },
                        [](py::tuple t) {
                            return mesh_from_arrays(t[0].cast<RArr>(), t[1].cast<py::array_t<long long, py::array::c_style | py::array::forcecast>>(),
                                                    false, Vec3());
                        }));

    // --- Generatoren ------------------------------------------------------------------------------------------------------------
    m.def("make_icosphere", &make_icosphere, "n"_a, "radius"_a = 1.0, "Ikosaederkugel mit 20 n^2 Dreiecken");
    m.def("make_icosphere_graded", &make_icosphere_graded, "n"_a, "pole"_a, "lam"_a, "radius"_a = 1.0,
          "Ikosaederkugel, konform zum Pol verdichtet (Elementgroesse am Pol etwa lam h)");
    m.def("make_cube_uniform", &make_cube_uniform, "n"_a, "Wuerfel [-1,1]^3, gleichmaessig, 12 n^2 Dreiecke");
    m.def("make_cube_graded", &make_cube_graded, "L"_a, "Wuerfel [-1,1]^3, zu den Kanten gradiert, 48 (L+1)^2 Dreiecke");
    m.def("translated", &translated, "mesh"_a, "shift"_a, "scale"_a = 1.0, "x -> scale x + shift");
    m.def("offset_surface", [](const TriangleMesh& me, real d) { return nogil([&] { return offset_surface(me, d); }); }, "mesh"_a, "d"_a,
          "Parallelflaeche im Abstand d (d > 0 nach aussen); RuntimeError bei Faltung oder Durchdringung");
    m.def("distance_to_surface", [](const TriangleMesh& me, const RArr& q, real rmax) {
              auto pts = to_points(q);
              auto d = nogil([&] { return distance_to_surface(me, pts, rmax); });
              return to_numpy(d);
          }, "mesh"_a, "points"_a, "rmax"_a, "kleinster Abstand der Punkte von der Flaeche (hoechstens rmax)");
    m.def("winding_number", [](const TriangleMesh& me, const RArr& x) -> py::object {
              auto pts = to_points(x);
              std::vector<real> w(pts.size());
              for (std::size_t i = 0; i < pts.size(); ++i) w[i] = winding_number(me, pts[i]);
              if (single_point(x)) return py::float_(w[0]);
              return to_numpy(w);
          }, "mesh"_a, "x"_a, "Windungszahl: 1 innen, 0 aussen (Punkt oder Punkte (M, 3))");
    m.def("signed_volume", &signed_volume, "mesh"_a, "Volumen mit Vorzeichen (> 0: nach aussen orientiert)");
    m.def("require_separated", &require_separated, "a"_a, "b"_a, "extra"_a = 0.0, "what"_a = std::string("Koerper"),
          "RuntimeError, wenn sich die Flaechen beruehren oder durchdringen");

    // --- Gmsh ---------------------------------------------------------------------------------------------------------------------
    m.def("read_gmsh", [](const std::string& path, real scale) {
              auto bodies = read_gmsh(path, scale);
              py::list out;
              for (auto& b : bodies) out.append(py::make_tuple(b.tag, std::move(b.mesh)));
              return out;
          }, "path"_a, "scale"_a = 1.0, "Gmsh 2.2/4.1 (ASCII): Liste von (tag, TriangleMesh) je Koerper, nach aussen orientiert");
    m.def("write_gmsh22", &write_gmsh22, "path"_a, "bodies"_a, "Gmsh 2.2 mit physikalischen Gruppen 1, 2, ...");

    // --- Richtungen ---------------------------------------------------------------------------------------------------------------
    m.def("lebedev", [](int n) {
              auto dirs = lebedev(n);
              std::vector<Vec3> d; std::vector<real> w;
              for (auto& x : dirs) { d.push_back(x.d); w.push_back(x.w); }
              return py::make_tuple(points_to_numpy(d), to_numpy(w));
          }, "npoints"_a, "Lebedev-Richtungen (1, 6, 14, 26): (Richtungen (n, 3), Gewichte (n,), Summe 1)");

    // --- Gruppierung fuer Blockvorkonditionierer ---------------------------------------------------------------------------------
    py::class_<FeatureSet>(m, "FeatureSet", "Ecken und Kanten eines Koerpers (Kanten-/Eckbloecke)")
        .def(py::init<>())
        .def_static("cube", &FeatureSet::cube, "a"_a = 1.0, "center"_a = Vec3(0, 0, 0), "Wuerfel center + [-a, a]^3")
        .def_readwrite("vertices", &FeatureSet::vertices)
        .def_readwrite("edges", &FeatureSet::edges)
        .def("append", &FeatureSet::append, "other"_a);
    m.def("group_by_features", &group_by_features, "mesh"_a, "features"_a, "R"_a, "Dreiecksgruppen um Ecken und Kanten");
    m.def("group_by_clusters", &group_by_clusters, "mesh"_a, "max_size"_a, "Dreiecksgruppen aus Clusterbaumblaettern");
}

}  // namespace cbem::py_bind
