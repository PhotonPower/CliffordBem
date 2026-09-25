// Streuung an beliebigen Koerpern aus einer Gmsh-Datei (ASCII 2.2/4.1). Ein Medium je Koerper (physikalische
// Gruppe), in der Reihenfolge aufsteigender Tags; fehlende Medien werden mit dem letzten angegebenen aufgefuellt.
// Beispiel: scatter_mesh --mesh stab.msh --omega 0.5 --media "-11,1.2;2.25,0" --chi "0;0.1" --pol circ --dir 0,0,1
#include <cstdio>
#include <sstream>
#include <string>
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/fields.hpp"
using namespace cbem;
static std::vector<std::string> split(const std::string& s, char c) { std::vector<std::string> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, c)) v.push_back(t); return v; }
static cplx cval(const std::string& s) { auto c = s.find(','); return {std::stod(s.substr(0, c)), c == std::string::npos ? 0.0 : std::stod(s.substr(c + 1))}; }
int main(int argc, char** argv) {
    std::string precond = "point";   // point | cluster:G (Bloecke auf Clustern mit <= G Dreiecken) | hodlr:eps[:leaf] (hierarchische Faktorisierung)
    std::string path, media = "2.25,0", chis = "0", pol = "lin"; double om = 1.0, scale = 1.0, heps = 1e-4, tol = 1e-6; Vec3 d(0, 0, 1);
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--precond") precond = nxt();
        else if (o == "--mesh") path = nxt(); else if (o == "--omega") om = std::stod(nxt()); else if (o == "--media") media = nxt();
        else if (o == "--chi") chis = nxt(); else if (o == "--pol") pol = nxt(); else if (o == "--scale") scale = std::stod(nxt());
        else if (o == "--dir") { auto v = split(nxt(), ','); d = Vec3(std::stod(v[0]), std::stod(v[1]), std::stod(v[2])); d = d / norm(d); }
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt());
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    if (path.empty()) { std::printf("--mesh fehlt\n"); return 1; }
    auto bodies = read_gmsh(path, scale);
    auto ms = split(media, ';'), cs = split(chis, ';');
    std::vector<TriangleMesh> parts; std::vector<Medium> med;
    for (std::size_t b = 0; b < bodies.size(); ++b) {
        Medium md; md.eps = cval(ms[std::min(b, ms.size() - 1)]); md.chi = cval(cs[std::min(b, cs.size() - 1)]);
        parts.push_back(bodies[b].mesh); med.push_back(md);
        std::printf("Koerper %zu (Tag %d): %zu Dreiecke, Volumen %.4f, eps = %g%+gi, chi = %g%+gi\n", b, bodies[b].tag, bodies[b].mesh.size(),
                    signed_volume(bodies[b].mesh), md.eps.real(), md.eps.imag(), md.chi.real(), md.chi.imag());
    }
    HMatrixParams hp; hp.eps = heps; SolveOptions so; so.tol = tol;
    ScatteringProblem P(parts, med, om, {}, hp);
    if (precond.rfind("hodlr:", 0) == 0) { auto q = precond.substr(6); auto c = q.find(':'); HodlrParams hpar; hpar.eps = std::stod(q.substr(0, c)); if (c != std::string::npos) hpar.leaf = std::stoul(q.substr(c + 1)); P.use_hodlr_preconditioner(hpar); std::printf("HODLR: eps %.0e, max. Rang %zu, %.0f MB, Aufbau %.1f s\n", hpar.eps, P.hodlr()->max_rank(), P.hodlr()->bytes() / 1048576.0, P.hodlr()->seconds()); }
        else if (precond.rfind("cluster:", 0) == 0) {
        P.use_block_preconditioner(group_by_clusters(P.mesh(), std::stoul(precond.substr(8))));
        std::printf("Blockvorkonditionierung: groesster Block %zu Unbekannte, Aufbau %.1f s\n", P.block_preconditioner()->max_block(), P.block_preconditioner()->seconds());
    }
    if (pol == "circ") {
        auto rp = P.solve_plane_wave(d, circular_polarization(d, +1), so), rm = P.solve_plane_wave(d, circular_polarization(d, -1), so);
        std::printf("sigma_ext: + %.6f, - %.6f, CD %.6f  (GMRES %d/%d, H %.0f MB)\n", rp.sigma_ext, rm.sigma_ext, rp.sigma_ext - rm.sigma_ext, rp.iterations, rm.iterations, P.hmatrix_bytes() / 1048576.0);
    } else {
        Vec3 a = std::abs(d.z) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0); Vec3 v = cross(d, a); v = v / norm(v); Vec3 u = cross(v, d);
        auto r = P.solve_plane_wave(d, CVec3{u.x, u.y, u.z}, so);
        std::printf("sigma_ext = %.6f  (GMRES %d, H %.0f MB)\n", r.sigma_ext, r.iterations, P.hmatrix_bytes() / 1048576.0);
    }
    return 0;
}
