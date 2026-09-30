// Zerfallsraten eines elektrischen Dipols (Fluorophor) vor Nanostrukturen (v0.26): gesamt, strahlend, nicht strahlend,
// veraenderte Quantenausbeute; Abstandsreihen, Orientierung radial / tangential / gemittelt.
// Geometrie: --sphere n (Radius = --unit nm, Kugelnetz zum Fusspunkt konform verdichtet, --graded auto | lambda | 1; auto:
// Elementgroesse am Fusspunkt ~ d/8, siehe docs/results_dipole.md) oder
// --mesh datei.msh (Dipolorte mit --pos). Schichten: --coating "d:Material;..." mit dem Zweitor (Dipol ausserhalb der Schichten).
// Abstaende: --dist "1,2,5,10" (nm vor der aeussersten Flaeche, entlang +x) oder --pos "x,y,z" (nm).
// Quantenausbeute: --q0 (intrinsisch, Standard 1): q = gamma_rad / (gamma_tot + (1 - q0)/q0), Raten relativ zum freien Dipol.
// Beispiel: dipole --sphere 12 --unit 20 --materials Au --nbg 1.33 --lambda 600 --dist "2,5,10,20" --orient average --q0 0.5
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "cbem/core/materials.hpp"
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/twoport_layer_problem.hpp"
#include "cbem/sources/dipole.hpp"
using namespace cbem;
static std::vector<std::string> split(const std::string& s, char c) { std::vector<std::string> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, c)) v.push_back(t); return v; }
int main(int argc, char** argv) {
    std::string mesh, mats = "Au", coating, dists, pos, orient = "radial", graded = "auto", csv, datadir = "data/materials";
    int sph = 0; double unit = 20, nbg = 1.0, lambda = 600, q0 = 1.0, heps = 1e-6, tol = 1e-9;
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--sphere") sph = std::stoi(nxt()); else if (o == "--mesh") mesh = nxt(); else if (o == "--unit") unit = std::stod(nxt());
        else if (o == "--materials") mats = nxt(); else if (o == "--nbg") nbg = std::stod(nxt()); else if (o == "--lambda") lambda = std::stod(nxt());
        else if (o == "--coating") coating = nxt(); else if (o == "--dist") dists = nxt(); else if (o == "--pos") pos = nxt();
        else if (o == "--orient") orient = nxt(); else if (o == "--graded") graded = nxt(); else if (o == "--q0") q0 = std::stod(nxt());
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt()); else if (o == "--csv") csv = nxt();
        else if (o == "--data") datadir = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    const real om = 2 * pi * unit / lambda;
    const Medium core{make_material(mats, datadir)->eps(lambda), 1.0, 0.0}, bg{nbg * nbg, 1.0, 0.0};
    std::vector<Coating> cs; real T = 0;
    for (auto& c : split(coating, ';')) { if (c.empty()) continue; auto q = split(c, ':');
        cs.push_back(Coating{std::stod(q[0]) / unit, Medium{make_material(q[1], datadir)->eps(lambda), 1.0, 0.0}}); T += cs.back().thickness; }
    // Dipolorte (Einheiten des Netzes)
    std::vector<Vec3> sites; std::vector<real> dnm;
    if (!pos.empty()) { auto v = split(pos, ','); sites.push_back(Vec3(std::stod(v[0]), std::stod(v[1]), std::stod(v[2])) / unit); dnm.push_back(-1); }
    for (auto& s : split(dists, ',')) { if (s.empty()) continue; const real d = std::stod(s); dnm.push_back(d); sites.push_back(Vec3(1.0 + T + d / unit, 0, 0)); }
    std::vector<std::pair<std::string, CVec3>> dirs;
    if (orient == "radial" || orient == "average") dirs.push_back({"radial", CVec3{1.0, 0.0, 0.0}});
    if (orient == "tangential" || orient == "average") dirs.push_back({"tangential", CVec3{0.0, 0.0, 1.0}});
    std::ofstream f; if (!csv.empty()) { f.open(csv, std::ios::app); f.seekp(0, std::ios::end); if (f.tellp() == 0) f << "lambda_nm,dist_nm,orientation,gamma_tot,gamma_rad,gamma_nr,q,iterations,lambda_graded,hmin_nm\n"; }   // anhaengen (Abstaende in getrennten Laeufen)
    std::printf("lambda %.1f nm, eps Kern %.3f%+.3fi, q0 = %.2f\n  d (nm)  Orientierung   gamma_tot   gamma_rad   gamma_nr        q    It.  (verdichtet, h_min nm)\n", lambda, core.eps.real(), core.eps.imag(), q0);
    for (std::size_t s = 0; s < sites.size(); ++s) {
        // Netz: Kugel zum Fusspunkt verdichtet (h am Fusspunkt ~ d/3), sonst Gmsh
        TriangleMesh m; real lam = 1;
        if (sph > 0) {
            // Faustregel aus der Validierung: Elementgroesse am Fusspunkt ~ d/8 fuer wenige Prozent in der Gesamtrate
            // (d/h = 3,4: +21 / +34 %, d/h = 8,6: -0,6 / +5,8 % radial / tangential); lambda nicht unter 0,2, sonst wird die
            // Gegenseite zu grob (strahlende Rate) -- dann n erhoehen
            const real d = std::max(1e-6, norm(sites[s]) - 1.0 - T), hu = 1.05 / sph, want = (d / 8) / hu;
            lam = graded == "auto" ? std::min(1.0, std::max(0.2, want)) : std::stod(graded);
            if (graded == "auto" && want < 0.2)
                std::printf("  Hinweis: fuer d = %.2f nm waere lambda = %.2f noetig; mit lambda = 0,2 ist h am Fusspunkt ~ d/%.1f -- n erhoehen\n",
                            d * unit, want, d / (0.2 * hu));
            m = lam < 1 ? make_icosphere_graded(sph, sites[s], lam) : make_icosphere(sph);
        } else m = read_gmsh(mesh)[0].mesh;
        real hmin = 1e9; for (real h : m.hmax) hmin = std::min(hmin, h);
        HMatrixParams hp; hp.eps = heps; SolveOptions so; so.tol = tol;
        std::unique_ptr<ScatteringProblem> PS; std::unique_ptr<TwoPortLayerProblem> P2; TriangleMesh outer;
        if (cs.empty()) { PS = std::make_unique<ScatteringProblem>(std::vector<TriangleMesh>{m}, std::vector<Medium>{core}, om, bg, hp); outer = PS->mesh(); }
        else { P2 = std::make_unique<TwoPortLayerProblem>(std::vector<TwoPortBody>{TwoPortBody{TwoPortLayerProblem::layer_surfaces(m, cs), core, cs}}, om, bg, hp); outer = P2->outer_mesh(); }
        if (winding_number(outer, sites[s]) > 0.5) { std::printf("Dipolort innerhalb des Koerpers\n"); return 1; }
        real avg_t = 0, avg_r = 0;
        for (auto& [name, p] : dirs) {
            const auto b = project_dipole(outer, bg, om, sites[s], p);
            std::vector<cplx> h; int its;
            if (PS) { auto r = PS->solve_rhs(b, so); h = r.h; its = r.iterations; }
            else { auto r = P2->solve_rhs(b, so); h.assign(r.h.begin(), r.h.begin() + b.size()); its = r.iterations; }
            const DipoleRates R = dipole_rates(outer, h, b, bg, om, sites[s], p);
            const real q = R.radiative / (R.total + (1 - q0) / q0);
            const real w = dirs.size() == 2 ? (name == "radial" ? 1.0 / 3 : 2.0 / 3) : 1.0;   // Kugel: 1 radiale, 2 tangentiale Richtungen
            avg_t += w * R.total; avg_r += w * R.radiative;
            std::printf("  %6.2f  %-12s %11.4f %11.4f %10.4f %8.4f %5d  (%.2f, %.2f)%s\n", dnm[s], name.c_str(), R.total, R.radiative, R.nonradiative, q, its, lam, hmin * unit, R.too_close ? "  zu nah" : "");
            if (f) f << lambda << ',' << dnm[s] << ',' << name << ',' << R.total << ',' << R.radiative << ',' << R.nonradiative << ',' << q << ',' << its << ',' << lam << ',' << hmin * unit << std::endl;
            std::fflush(stdout);
        }
        if (dirs.size() == 2) {
            const real q = avg_r / (avg_t + (1 - q0) / q0);
            std::printf("  %6.2f  %-12s %11.4f %11.4f %10.4f %8.4f\n", dnm[s], "gemittelt", avg_t, avg_r, avg_t - avg_r, q);
            if (f) f << lambda << ',' << dnm[s] << ",average," << avg_t << ',' << avg_r << ',' << avg_t - avg_r << ',' << q << ",0," << lam << ',' << hmin * unit << std::endl;
            std::fflush(stdout);
        }
    }
    return 0;
}
