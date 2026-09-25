// Spektren: Extinktion (und Zirkulardichroismus) ueber der Vakuumwellenlaenge, mit dispersiven Materialien
// (Johnson-Christy Au/Ag) in einem Hintergrundmedium, optional orientierungsgemittelt (Lebedev 6/14/26).
// Laengeneinheit der Geometrie: --unit nm (Kugel: Radius 1 Einheit; Gmsh-Datei: Koordinaten in Einheiten).
// Beispiele:
//   spectrum --sphere 12 --unit 40 --materials Au --nbg 1.33 --lambda 450:650:10 --csv results/au40_water.csv
//   spectrum --mesh examples/bornkuhn_60.msh --unit 20 --materials Au --lambda 500:900:25 --pol circ --orient 14
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "cbem/core/materials.hpp"
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/fields.hpp"
#include "cbem/sources/orientation.hpp"
using namespace cbem;
static std::vector<std::string> split(const std::string& s, char c) { std::vector<std::string> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, c)) v.push_back(t); return v; }
int main(int argc, char** argv) {
    std::string precond = "point";   // point | cluster:G (Bloecke auf Clustern mit <= G Dreiecken) | hodlr:eps[:leaf] (hierarchische Faktorisierung)
    std::string mesh, mats = "Au", chis = "0", pol = "lin", lam = "500:600:50", csv, datadir = "data/materials";
    int sph = 0, orient = 1; bool verbose = false; double unit = 1.0, nbg = 1.0, heps = 1e-4, tol = 1e-6;
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--precond") precond = nxt();
        else if (o == "--mesh") mesh = nxt(); else if (o == "--sphere") sph = std::stoi(nxt()); else if (o == "--unit") unit = std::stod(nxt());
        else if (o == "--materials") mats = nxt(); else if (o == "--chi") chis = nxt(); else if (o == "--nbg") nbg = std::stod(nxt());
        else if (o == "--lambda") lam = nxt(); else if (o == "--pol") pol = nxt(); else if (o == "--orient") orient = std::stoi(nxt());
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt());
        else if (o == "--csv") csv = nxt(); else if (o == "--data") datadir = nxt(); else if (o == "--verbose") verbose = true;
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    std::vector<TriangleMesh> parts;
    if (sph > 0) parts.push_back(make_icosphere(sph));
    else { for (auto& b : read_gmsh(mesh)) parts.push_back(b.mesh); }
    auto ms = split(mats, ';'), cs = split(chis, ';');
    std::vector<std::shared_ptr<Material>> mat; std::vector<cplx> chi;
    for (std::size_t b = 0; b < parts.size(); ++b) {
        mat.push_back(make_material(ms[std::min(b, ms.size() - 1)], datadir));
        auto c = cs[std::min(b, cs.size() - 1)]; auto q = c.find(',');
        chi.push_back(cplx(std::stod(c.substr(0, q)), q == std::string::npos ? 0.0 : std::stod(c.substr(q + 1))));
    }
    std::vector<double> lams; { auto r = split(lam, ':'); if (r.size() == 3) for (double l = std::stod(r[0]); l <= std::stod(r[1]) + 1e-9; l += std::stod(r[2])) lams.push_back(l); else for (auto& t : split(lam, ',')) lams.push_back(std::stod(t)); }
    const auto dirs = lebedev(orient);
    std::ofstream f; if (!csv.empty()) { f.open(csv, std::ios::app); f.seekp(0, std::ios::end); if (f.tellp() == 0) f << "lambda_nm,unit_nm,nbg,N,orient,pol,sigma_nm2,sigma_plus_nm2,sigma_minus_nm2,CD_nm2,iterations,t_s\n"; }
    std::printf("%s: %zu Koerper, %zu Dreiecke gesamt; Einheit %.3g nm, n_Hintergrund %.3f, %zu Richtung(en)\n",
                sph ? "Kugel" : mesh.c_str(), parts.size(), [&] { std::size_t n = 0; for (auto& p : parts) n += p.size(); return n; }(), unit, nbg, dirs.size());
    std::printf("%9s %13s %13s %13s %13s %6s %7s\n", "lambda nm", "sigma nm^2", "sigma+ nm^2", "sigma- nm^2", "CD nm^2", "It.", "Zeit s");
    for (double L : lams) {
        auto t0 = std::chrono::steady_clock::now();
        const real om = 2 * pi * unit / L;                     // k0 in 1/Einheit
        std::vector<Medium> med;
        for (std::size_t b = 0; b < parts.size(); ++b) med.push_back(Medium{mat[b]->eps(L), 1.0, chi[b]});
        HMatrixParams hp; hp.eps = heps; SolveOptions so; so.tol = tol;
        ScatteringProblem P(parts, med, om, Medium{nbg * nbg, 1.0, 0.0}, hp);
        if (precond.rfind("hodlr:", 0) == 0) { auto q = precond.substr(6); auto c = q.find(':'); HodlrParams hpar; hpar.eps = std::stod(q.substr(0, c)); if (c != std::string::npos) hpar.leaf = std::stoul(q.substr(c + 1)); P.use_hodlr_preconditioner(hpar); std::printf("HODLR: eps %.0e, max. Rang %zu, %.0f MB, Aufbau %.1f s\n", hpar.eps, P.hodlr()->max_rank(), P.hodlr()->bytes() / 1048576.0, P.hodlr()->seconds()); }
        else if (precond.rfind("cluster:", 0) == 0) P.use_block_preconditioner(group_by_clusters(P.mesh(), std::stoul(precond.substr(8))));
        real s = 0, sp = 0, sm = 0; int its = 0;
        for (const auto& dw : dirs) {
            if (pol == "circ") {
                auto a = P.solve_plane_wave(dw.d, circular_polarization(dw.d, +1), so), b = P.solve_plane_wave(dw.d, circular_polarization(dw.d, -1), so);
                sp += dw.w * a.sigma_ext; sm += dw.w * b.sigma_ext; its = std::max(its, std::max(a.iterations, b.iterations));
                if (verbose) std::printf("      d = (%+.3f, %+.3f, %+.3f): sigma+ %.5g, sigma- %.5g, CD %+.4g nm^2\n", dw.d.x, dw.d.y, dw.d.z,
                                         a.sigma_ext * unit * unit, b.sigma_ext * unit * unit, (a.sigma_ext - b.sigma_ext) * unit * unit);
            } else {
                // linear: Mittel ueber zwei orthogonale Polarisationen = unpolarisiert (fuer orient > 1), sonst p = u
                CVec3 u = circular_polarization(dw.d, +1), v = circular_polarization(dw.d, -1);
                CVec3 pu{0.5 * (u[0] + v[0]), 0.5 * (u[1] + v[1]), 0.5 * (u[2] + v[2])};
                auto a = P.solve_plane_wave(dw.d, pu, so); real sa = a.sigma_ext; its = std::max(its, a.iterations);
                if (dirs.size() > 1) {
                    CVec3 pv{(u[0] - v[0]) / cplx(0, 2), (u[1] - v[1]) / cplx(0, 2), (u[2] - v[2]) / cplx(0, 2)};
                    sa = 0.5 * (sa + P.solve_plane_wave(dw.d, pv, so).sigma_ext);
                }
                s += dw.w * sa;
            }
        }
        if (pol == "circ") s = 0.5 * (sp + sm);
        const real u2 = unit * unit;
        double ts = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("%9.1f %13.5g %13.5g %13.5g %13.5g %6d %7.1f\n", L, s * u2, sp * u2, sm * u2, (sp - sm) * u2, its, ts); std::fflush(stdout);
        if (f) { f << L << ',' << unit << ',' << nbg << ',' << P.mesh().size() << ',' << dirs.size() << ',' << pol << ',' << s * u2 << ',' << sp * u2 << ',' << sm * u2 << ','
                   << (sp - sm) * u2 << ',' << its << ',' << ts << '\n'; f.flush(); }
    }
    return 0;
}
