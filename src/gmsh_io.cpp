#include "cbem/geometry/gmsh_io.hpp"
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace cbem {

real signed_volume(const TriangleMesh& m) {
    real v = 0;
    for (const auto& t : m.T) v += dot(m.P[t[0]], cross(m.P[t[1]], m.P[t[2]]));
    return v / 6.0;
}

namespace {
// Knotenzahl je Gmsh-Elementtyp (zum Ueberspringen); 0: unbekannt
int gmsh_nodes(int type) {
    switch (type) {
        case 1: return 2;  case 2: return 3;  case 3: return 4;  case 4: return 4;  case 5: return 8;  case 6: return 6;
        case 7: return 5;  case 8: return 3;  case 9: return 6;  case 10: return 9; case 11: return 10; case 12: return 27;
        case 13: return 18; case 14: return 14; case 15: return 1; case 16: return 8; case 17: return 20; case 18: return 15;
        case 19: return 13; case 20: return 9; case 21: return 10; case 26: return 4; case 27: return 5; case 29: return 20;
        default: return 0;
    }
}

// Knoten global einlesen, Dreiecke je Koerper sammeln (erste oder zweite Ordnung), dann kompakt umnummerieren
struct Builder {
    std::map<long, Vec3> nodes;
    std::map<int, std::vector<std::array<long, 6>>> tris;   // Ecken 0..2, Kantenmitten 3..5 (-1: erste Ordnung)
    void add(int body, int type, const long* v) {
        std::array<long, 6> t{v[0], v[1], v[2], -1, -1, -1};
        if (type == 9) { t[3] = v[3]; t[4] = v[4]; t[5] = v[5]; }
        tris[body].push_back(t);
    }
    const Vec3& node(long id) const {
        auto nd = nodes.find(id); if (nd == nodes.end()) throw std::runtime_error("Gmsh: unbekannter Knoten");
        return nd->second;
    }
    // ebene Netze (nur Ecken)
    std::vector<GmshBody> finish(real scale) const {
        std::vector<GmshBody> out;
        for (auto& [tag, list] : tris) {
            GmshBody b; b.tag = tag; std::map<long, int> id;
            for (auto& t : list) {
                std::array<int, 3> loc{};
                for (int q = 0; q < 3; ++q) {
                    auto it = id.find(t[q]);
                    if (it == id.end()) {
                        int k = static_cast<int>(b.mesh.P.size()); b.mesh.P.push_back(node(t[q]) * scale); id[t[q]] = k; loc[q] = k;
                    } else loc[q] = it->second;
                }
                b.mesh.T.push_back(loc);
            }
            b.mesh.compute_geometry();
            if (signed_volume(b.mesh) < 0) { for (auto& t : b.mesh.T) std::swap(t[1], t[2]); b.mesh.compute_geometry(); }
            out.push_back(std::move(b));
        }
        return out;
    }
    // quadratische Netze: Kantenmitten aus den Knoten 3..5; beim Umorientieren (v0, v2, v1) werden die Mitten zu (M20, M12, M01)
    std::vector<GmshQuadraticBody> finish_quadratic(real scale) const {
        std::vector<GmshQuadraticBody> out;
        for (auto& [tag, list] : tris) {
            GmshQuadraticBody b; b.tag = tag; std::map<long, int> id;
            TriangleMesh& f = b.mesh.flat;
            for (auto& t : list) {
                if (t[3] < 0) throw std::runtime_error("Gmsh: Koerper " + std::to_string(tag) + " enthaelt Dreiecke erster Ordnung "
                                                       "(mit Mesh.ElementOrder = 2 vernetzen oder make_quadratic verwenden)");
                std::array<int, 3> loc{};
                for (int q = 0; q < 3; ++q) {
                    auto it = id.find(t[q]);
                    if (it == id.end()) { int k = static_cast<int>(f.P.size()); f.P.push_back(node(t[q]) * scale); id[t[q]] = k; loc[q] = k; }
                    else loc[q] = it->second;
                }
                f.T.push_back(loc);
                b.mesh.mid.push_back({node(t[3]) * scale, node(t[4]) * scale, node(t[5]) * scale});
            }
            f.compute_geometry();
            if (signed_volume(f) < 0) {
                for (std::size_t e = 0; e < f.T.size(); ++e) { std::swap(f.T[e][1], f.T[e][2]); auto& M = b.mesh.mid[e]; M = {M[2], M[1], M[0]}; }
                f.compute_geometry();
            }
            b.mesh.compute_geometry();
            out.push_back(std::move(b));
        }
        return out;
    }
};

Builder parse_gmsh(const std::string& path) {
    std::ifstream in(path); if (!in) throw std::runtime_error("Gmsh: Datei nicht lesbar: " + path);
    Builder B; std::string line; double version = 0;
    std::map<int, int> entity_phys;                                      // 4.1: Flaechen-Entitaet -> physikalischer Tag
    while (std::getline(in, line)) {
        if (line == "$MeshFormat") { std::getline(in, line); std::istringstream s(line); int ft; s >> version >> ft; if (ft != 0) throw std::runtime_error("Gmsh: nur ASCII"); }
        else if (line == "$Entities" && version >= 4) {
            std::getline(in, line); std::istringstream s(line); long np, nc, ns, nv; s >> np >> nc >> ns >> nv;
            for (long i = 0; i < np + nc; ++i) std::getline(in, line);
            for (long i = 0; i < ns; ++i) {
                std::getline(in, line); std::istringstream e(line); int tag; double x; e >> tag; for (int q = 0; q < 6; ++q) e >> x;
                int nph; e >> nph; int ph = tag; if (nph > 0) e >> ph; entity_phys[tag] = ph;
            }
        }
        else if (line == "$Nodes") {
            std::getline(in, line); std::istringstream s(line);
            if (version < 4) { long n; s >> n; for (long i = 0; i < n; ++i) { long id; double x, y, z; in >> id >> x >> y >> z; B.nodes[id] = Vec3(x, y, z); } }
            else {
                long nb, nn, mn, mx; s >> nb >> nn >> mn >> mx;
                for (long b = 0; b < nb; ++b) {
                    int dim, tag, par; long cnt; in >> dim >> tag >> par >> cnt;
                    // parametrische Knoten (Mesh.SaveParametric): u auf Kurven, u v auf Flaechen; werden uebersprungen
                    const int npar = (par != 0 && (dim == 1 || dim == 2)) ? dim : 0;
                    std::vector<long> ids(cnt); for (auto& v : ids) in >> v;
                    for (long i = 0; i < cnt; ++i) {
                        double x, y, z, u; in >> x >> y >> z; for (int q = 0; q < npar; ++q) in >> u;
                        B.nodes[ids[i]] = Vec3(x, y, z);
                    }
                }
            }
        }
        else if (line == "$Elements") {
            std::getline(in, line); std::istringstream s(line);
            if (version < 4) {
                long n; s >> n; std::getline(in, line);
                for (long i = 0; i < n; ++i) {
                    if (i > 0) std::getline(in, line);
                    std::istringstream e(line); long id; int type, ntag; e >> id >> type >> ntag;
                    std::vector<int> tags(ntag); for (auto& t : tags) e >> t;
                    if (type != 2 && type != 9) continue;
                    long v[6] = {0, 0, 0, 0, 0, 0}; for (int q = 0; q < (type == 9 ? 6 : 3); ++q) e >> v[q];
                    B.add(ntag > 0 ? tags[0] : 0, type, v);
                }
            } else {
                long nb, ne, mn, mx; s >> nb >> ne >> mn >> mx;
                for (long b = 0; b < nb; ++b) {
                    int dim, tag, type; long cnt; in >> dim >> tag >> type >> cnt;
                    const int nper = gmsh_nodes(type);
                    if (nper == 0) throw std::runtime_error("Gmsh: nicht unterstuetzter Elementtyp " + std::to_string(type));
                    for (long i = 0; i < cnt; ++i) {
                        long id; in >> id; long v[27] = {};
                        for (int q = 0; q < nper; ++q) in >> v[q];
                        if ((type == 2 || type == 9) && dim == 2) B.add(entity_phys.count(tag) ? entity_phys[tag] : tag, type, v);
                    }
                }
            }
        }
    }
    if (B.tris.empty()) throw std::runtime_error("Gmsh: keine Dreiecke gefunden");
    return B;
}
}  // namespace

std::vector<GmshBody> read_gmsh(const std::string& path, real scale) { return parse_gmsh(path).finish(scale); }

std::vector<GmshQuadraticBody> read_gmsh_quadratic(const std::string& path, real scale) { return parse_gmsh(path).finish_quadratic(scale); }

void write_gmsh22(const std::string& path, const std::vector<TriangleMesh>& bodies) {
    std::ofstream o(path); o.precision(17);
    o << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
    std::size_t nn = 0, ne = 0; for (auto& b : bodies) { nn += b.P.size(); ne += b.T.size(); }
    o << nn << "\n"; std::size_t off = 0;
    for (auto& b : bodies) { for (std::size_t i = 0; i < b.P.size(); ++i) o << off + i + 1 << ' ' << b.P[i].x << ' ' << b.P[i].y << ' ' << b.P[i].z << "\n"; off += b.P.size(); }
    o << "$EndNodes\n$Elements\n" << ne << "\n"; off = 0; std::size_t id = 1;
    for (std::size_t k = 0; k < bodies.size(); ++k) {
        for (auto& t : bodies[k].T) o << id++ << " 2 2 " << k + 1 << ' ' << k + 1 << ' ' << off + t[0] + 1 << ' ' << off + t[1] + 1 << ' ' << off + t[2] + 1 << "\n";
        off += bodies[k].P.size();
    }
    o << "$EndElements\n";
}

void write_gmsh22_quadratic(const std::string& path, const std::vector<QuadraticMesh>& bodies) {
    std::ofstream o(path); o.precision(17);
    // Knoten: Ecken je Koerper, dann je Kante eine Mitte (Schluessel: sortiertes Eckenpaar)
    std::vector<Vec3> P; std::vector<std::vector<std::array<std::size_t, 6>>> el(bodies.size());
    for (std::size_t k = 0; k < bodies.size(); ++k) {
        const QuadraticMesh& q = bodies[k];
        const std::size_t off = P.size();
        P.insert(P.end(), q.flat.P.begin(), q.flat.P.end());
        std::map<std::pair<int, int>, std::size_t> edge;
        for (std::size_t t = 0; t < q.size(); ++t) {
            std::array<std::size_t, 6> e{};
            for (int c = 0; c < 3; ++c) e[c] = off + q.flat.T[t][c];
            for (int c = 0; c < 3; ++c) {
                const int a = q.flat.T[t][c], b = q.flat.T[t][(c + 1) % 3];
                const auto key = std::make_pair(std::min(a, b), std::max(a, b));
                auto it = edge.find(key);
                if (it == edge.end()) { edge[key] = P.size(); e[3 + c] = P.size(); P.push_back(q.mid[t][c]); }
                else e[3 + c] = it->second;
            }
            el[k].push_back(e);
        }
    }
    o << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n" << P.size() << "\n";
    for (std::size_t i = 0; i < P.size(); ++i) o << i + 1 << ' ' << P[i].x << ' ' << P[i].y << ' ' << P[i].z << "\n";
    std::size_t ne = 0; for (auto& v : el) ne += v.size();
    o << "$EndNodes\n$Elements\n" << ne << "\n";
    std::size_t id = 1;
    for (std::size_t k = 0; k < el.size(); ++k)
        for (auto& e : el[k]) {
            o << id++ << " 9 2 " << k + 1 << ' ' << k + 1;
            for (int c = 0; c < 6; ++c) o << ' ' << e[c] + 1;
            o << "\n";
        }
    o << "$EndElements\n";
}

}  // namespace cbem
