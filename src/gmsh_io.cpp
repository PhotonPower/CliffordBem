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
// Knoten global einlesen, Dreiecke je Koerper sammeln, dann kompakt umnummerieren
struct Builder {
    std::map<long, Vec3> nodes;
    std::map<int, std::vector<std::array<long, 3>>> tris;
    std::vector<GmshBody> finish(real scale) {
        std::vector<GmshBody> out;
        for (auto& [tag, list] : tris) {
            GmshBody b; b.tag = tag; std::map<long, int> id;
            for (auto& t : list) {
                std::array<int, 3> loc{};
                for (int q = 0; q < 3; ++q) {
                    auto it = id.find(t[q]);
                    if (it == id.end()) {
                        auto nd = nodes.find(t[q]); if (nd == nodes.end()) throw std::runtime_error("Gmsh: unbekannter Knoten");
                        int k = static_cast<int>(b.mesh.P.size()); b.mesh.P.push_back(nd->second * scale); id[t[q]] = k; loc[q] = k;
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
};
}  // namespace

std::vector<GmshBody> read_gmsh(const std::string& path, real scale) {
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
                    std::vector<long> ids(cnt); for (auto& v : ids) in >> v;
                    for (long i = 0; i < cnt; ++i) { double x, y, z; in >> x >> y >> z; B.nodes[ids[i]] = Vec3(x, y, z); }
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
                    if (type != 2) continue;
                    long a, b, c; e >> a >> b >> c; B.tris[ntag > 0 ? tags[0] : 0].push_back({a, b, c});
                }
            } else {
                long nb, ne, mn, mx; s >> nb >> ne >> mn >> mx;
                for (long b = 0; b < nb; ++b) {
                    int dim, tag, type; long cnt; in >> dim >> tag >> type >> cnt;
                    int nper = type == 2 ? 3 : (type == 1 ? 2 : (type == 15 ? 1 : (type == 4 ? 4 : (type == 3 ? 4 : 0))));
                    if (nper == 0 && type != 2) throw std::runtime_error("Gmsh: nicht unterstuetzter Elementtyp");
                    for (long i = 0; i < cnt; ++i) {
                        long id; in >> id; std::array<long, 4> v{};
                        for (int q = 0; q < nper; ++q) in >> v[q];
                        if (type == 2 && dim == 2) { int ph = entity_phys.count(tag) ? entity_phys[tag] : tag; B.tris[ph].push_back({v[0], v[1], v[2]}); }
                    }
                }
            }
        }
    }
    if (B.tris.empty()) throw std::runtime_error("Gmsh: keine Dreiecke gefunden");
    return B.finish(scale);
}

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

}  // namespace cbem
