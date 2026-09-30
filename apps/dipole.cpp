// Zerfallsraten und Fluoreszenzverstaerkung eines Emitters (Fluorophor) vor Nanostrukturen (v0.26/v0.27).
// Geometrie: --sphere n (Radius = --unit nm) [--sphere-dimer g: zwei Kugeln entlang x, Spalt g nm zwischen den Kernen] oder
// --mesh datei.msh; Kugelnetze zum Fusspunkt des Emitters konform verdichtet (--graded auto | lambda | 1; auto: h ~ d/8).
// Schichten: --coating "d:Material;..." mit dem Zweitor (Emitter ausserhalb der Schichten).
// Emitterorte: --dist "1,2,5" (nm vor der aeussersten Flaeche, entlang +x; Kugel) oder --pos "x,y,z" (nm); Dimer: Spaltmitte.
// Emission bei --lambda (nm): Raten fuer die Orientierungen x, y, z (Kugel und Dimer mit Emitter auf der x-Achse: y = z).
// Anregung (optional) bei --lambda-exc (nm): ebene Welle --exc-dir (Standard z) mit linearer Polarisation --exc-pol (Standard x),
// lokales Feld am Emitterort; Fluoreszenzverstaerkung fuer fest, aber zufaellig orientierte Emitter:
//   F/F0 = sum_a |E_a|^2 q_a / (|E0|^2 q0),   q_a = gamma_rad,a / (gamma_tot,a + (1 - q0)/q0),
// Anregung und Quantenausbeute gemeinsam gemittelt (nicht das Produkt der Mittelwerte). --q0: intrinsische Quantenausbeute.
// Chiraler Emitter (v0.32): --chiral kappa -- zu jedem elektrischen Dipol p ein magnetischer m = i kappa sqrt(mu/eps) p
// (parallele Uebergangsdipole); ausgegeben werden g_lum = 2 (P+ - P-)/(P+ + P-) mit und ohne Nanostruktur sowie das Mittel ueber
// die Orientierungen aus den gemittelten Leistungen je Helizitaet (Beschriftung s wie circular_polarization und der CD).
// Beispiel: dipole --sphere 12 --unit 20 --materials Au --nbg 1.33 --lambda 600 --lambda-exc 580 --dist "2,5,10,20" --q0 0.1
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include "cbem/core/materials.hpp"
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/twoport_layer_problem.hpp"
#include "cbem/sources/dipole.hpp"
using namespace cbem;
static std::vector<std::string> split(const std::string& s, char c) { std::vector<std::string> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, c)) v.push_back(t); return v; }
static Vec3 vec3(const std::string& s) { auto v = split(s, ','); return Vec3(std::stod(v[0]), std::stod(v[1]), std::stod(v[2])); }

struct Geometry { std::vector<TriangleMesh> bodies; real lam = 1, hmin = 0; };

int main(int argc, char** argv) {
    std::string mesh, mats = "Au", coating, dists, pos, graded = "auto", csv, datadir = "data/materials", edir = "0,0,1", epol = "1,0,0";
    int sph = 0; double unit = 20, nbg = 1.0, lambda = 600, lexc = -1, q0 = 1.0, gap = -1, heps = 1e-6, tol = 1e-9, kappa = 0;
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--sphere") sph = std::stoi(nxt()); else if (o == "--sphere-dimer") gap = std::stod(nxt()); else if (o == "--mesh") mesh = nxt();
        else if (o == "--unit") unit = std::stod(nxt()); else if (o == "--materials") mats = nxt(); else if (o == "--nbg") nbg = std::stod(nxt());
        else if (o == "--lambda") lambda = std::stod(nxt()); else if (o == "--lambda-exc") lexc = std::stod(nxt());
        else if (o == "--exc-dir") edir = nxt(); else if (o == "--exc-pol") epol = nxt();
        else if (o == "--coating") coating = nxt(); else if (o == "--dist") dists = nxt(); else if (o == "--pos") pos = nxt();
        else if (o == "--graded") graded = nxt(); else if (o == "--q0") q0 = std::stod(nxt()); else if (o == "--chiral") kappa = std::stod(nxt());
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt()); else if (o == "--csv") csv = nxt();
        else if (o == "--data") datadir = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    auto coats_at = [&](real L, real& T) {
        std::vector<Coating> cs; T = 0;
        for (auto& c : split(coating, ';')) { if (c.empty()) continue; auto q = split(c, ':');
            cs.push_back(Coating{std::stod(q[0]) / unit, Medium{make_material(q[1], datadir)->eps(L), 1.0, 0.0}}); T += cs.back().thickness; }
        return cs;
    };
    real T = 0; coats_at(lambda, T);
    const real sx = gap >= 0 ? 1.0 + 0.5 * gap / unit : 0.0;                   // Dimer: Kugelmitten bei -sx, +sx
    std::vector<Vec3> sites; std::vector<real> dnm;
    if (!pos.empty()) { sites.push_back(vec3(pos) / unit); dnm.push_back(-1); }
    else if (gap >= 0) { sites.push_back(Vec3(0, 0, 0)); dnm.push_back(0.5 * gap - T * unit); }
    for (auto& s : split(dists, ',')) { if (s.empty() || gap >= 0) continue; const real d = std::stod(s); dnm.push_back(d); sites.push_back(Vec3(1.0 + T + d / unit, 0, 0)); }
    // Netz je Emitterort: Kugel(n) zum Fusspunkt verdichtet (h ~ d/8, lambda >= 0,2)
    auto geometry = [&](const Vec3& site) {
        Geometry G;
        if (sph > 0) {
            const real hu = 1.05 / sph;
            const std::vector<Vec3> centers = gap >= 0 ? std::vector<Vec3>{Vec3(-sx, 0, 0), Vec3(sx, 0, 0)} : std::vector<Vec3>{Vec3(0, 0, 0)};
            real d = 1e9; for (auto& c : centers) d = std::min(d, norm(site - c) - 1.0 - T);
            const real want = (std::max(d, 1e-6) / 8) / hu;
            G.lam = graded == "auto" ? std::min(1.0, std::max(0.2, want)) : std::stod(graded);
            if (graded == "auto" && want < 0.2) std::printf("  Hinweis: fuer d = %.2f nm waere lambda = %.2f noetig (h am Fusspunkt ~ d/%.1f) -- n erhoehen\n", d * unit, want, d / (0.2 * hu));
            for (auto& c : centers) G.bodies.push_back(translated(G.lam < 1 ? make_icosphere_graded(sph, site - c, G.lam) : make_icosphere(sph), c));
        } else for (auto& b : read_gmsh(mesh)) G.bodies.push_back(b.mesh);
        G.hmin = 1e9; for (auto& b : G.bodies) { TriangleMesh m = b; m.compute_geometry(); for (real h : m.hmax) G.hmin = std::min(G.hmin, h); }
        return G;
    };
    // Loeser fuer eine Wellenlaenge (Kern und Schichten dispersiv)
    struct Solver { std::unique_ptr<ScatteringProblem> PS; std::unique_ptr<TwoPortLayerProblem> P2; TriangleMesh outer; real om = 0; Medium bg; };
    auto make_solver = [&](const Geometry& G, real L) {
        Solver S; S.om = 2 * pi * unit / L; S.bg = Medium{nbg * nbg, 1.0, 0.0};
        const Medium core{make_material(mats, datadir)->eps(L), 1.0, 0.0}; real TT; const auto cs = coats_at(L, TT);
        HMatrixParams hp; hp.eps = heps;
        if (cs.empty()) { S.PS = std::make_unique<ScatteringProblem>(G.bodies, std::vector<Medium>(G.bodies.size(), core), S.om, S.bg, hp); S.outer = S.PS->mesh(); }
        else { std::vector<TwoPortBody> tb; for (auto& b : G.bodies) tb.push_back(TwoPortBody{TwoPortLayerProblem::layer_surfaces(b, cs), core, cs});
               S.P2 = std::make_unique<TwoPortLayerProblem>(tb, S.om, S.bg, hp); S.outer = S.P2->outer_mesh(); }
        return S;
    };
    auto solve = [&](const Solver& S, const std::vector<cplx>& b, int& its) {
        SolveOptions so; so.tol = tol; std::vector<cplx> h;
        if (S.PS) { auto r = S.PS->solve_rhs(b, so); h = r.h; its = r.iterations; }
        else { auto r = S.P2->solve_rhs(b, so); h.assign(r.h.begin(), r.h.begin() + b.size()); its = r.iterations; }
        return h;
    };
    std::ofstream f;
    if (!csv.empty()) { f.open(csv, std::ios::app); f.seekp(0, std::ios::end);
        if (f.tellp() == 0) f << "lambda_nm,lambda_exc_nm,dist_nm,orientation,gamma_tot,gamma_rad,gamma_nr,q,exc,iterations,lambda_graded,hmin_nm,kappa,rad_plus,rad_minus,glum,glum_free\n"; }
    std::printf("Emission %.1f nm", lambda); if (lexc > 0) std::printf(", Anregung %.1f nm", lexc); std::printf(", q0 = %.2f\n", q0);
    std::printf("  d (nm)  Orient.   gamma_tot   gamma_rad    gamma_nr        q   |E_a|^2/|E0|^2   It.  (lambda, h_min nm)\n");
    const char* NM[3] = {"x", "y", "z"};
    for (std::size_t s = 0; s < sites.size(); ++s) {
        const Geometry G = geometry(sites[s]);
        const bool axial = sph > 0 && std::abs(sites[s].y) < 1e-12 && std::abs(sites[s].z) < 1e-12;   // Emitter auf der x-Achse: y = z
        // Anregung: lokales Feld der ebenen Welle am Emitterort (|E0| = 1)
        real exc[3] = {-1, -1, -1};
        if (lexc > 0) {
            const Solver X = make_solver(G, lexc); const Vec3 d = vec3(edir) / norm(vec3(edir)); const Vec3 pv = vec3(epol);
            const CVec3 p{pv.x, pv.y, pv.z};
            const auto b = project_plane_wave(X.outer, X.bg.k(X.om), X.bg.eps, d, p); int its;
            const auto h = solve(X, b, its);
            const auto nf = exterior_near_field(X.outer, h, X.bg, X.om, d, p, {sites[s]});
            const real p2 = std::norm(p[0]) + std::norm(p[1]) + std::norm(p[2]);
            for (int a = 0; a < 3; ++a) exc[a] = std::norm(nf[0].E[a]) / p2;
        }
        // Emission: Raten je Orientierung
        const Solver S = make_solver(G, lambda);
        if (winding_number(S.outer, sites[s]) > 0.5) { std::printf("Emitterort innerhalb eines Koerpers\n"); return 1; }
        DipoleRates R[3]; int IT[3] = {0, 0, 0};
        for (int a = 0; a < 3; ++a) {
            if (axial && a == 1) continue;                                        // y wie z
            CVec3 p{0.0, 0.0, 0.0}; p[a] = 1.0;
            CVec3 md{0.0, 0.0, 0.0}; md[a] = cplx(0, kappa * std::sqrt(std::real(S.bg.mu / S.bg.eps)));   // chiraler Emitter
            const auto b = project_dipole(S.outer, S.bg, S.om, sites[s], p, 24, md);
            const auto h = solve(S, b, IT[a]);
            R[a] = dipole_rates(S.outer, h, b, S.bg, S.om, sites[s], p, 40, md);
        }
        if (axial) { R[1] = R[2]; IT[1] = IT[2]; }
        real F = 0, gt = 0, gr = 0;
        for (int a = 0; a < 3; ++a) {
            const real q = R[a].radiative / (R[a].total + (1 - q0) / q0);
            gt += R[a].total / 3; gr += R[a].radiative / 3;
            (void)F;
            std::printf("  %6.2f  %-7s %11.4f %11.4f %11.4f %8.4f %14.3f %6d  (%.2f, %.2f)%s", dnm[s], NM[a], R[a].total, R[a].radiative, R[a].nonradiative, q,
                        exc[a], IT[a], G.lam, G.hmin * unit, R[a].too_close ? "  zu nah" : "");
            if (kappa != 0) std::printf("   g_lum %+.5f (frei %+.5f)", R[a].glum, R[a].glum_free);
            std::printf("\n");
            if (f) f << lambda << ',' << lexc << ',' << dnm[s] << ',' << NM[a] << ',' << R[a].total << ',' << R[a].radiative << ',' << R[a].nonradiative << ',' << q << ','
                     << exc[a] << ',' << IT[a] << ',' << G.lam << ',' << G.hmin * unit << ',' << kappa << ',' << R[a].rad_plus << ',' << R[a].rad_minus
                     << ',' << R[a].glum << ',' << R[a].glum_free << std::endl;
        }
        if (lexc > 0) F = fluorescence_enhancement(exc, R, q0);                   // sum_a |E_a|^2 q_a / (|E0|^2 q0)
        std::printf("  %6.2f  Mittel  %11.4f %11.4f %11.4f %8.4f", dnm[s], gt, gr, gt - gr, gr / (gt + (1 - q0) / q0));
        if (lexc > 0) std::printf("   F/F0 = %.3f (fest, zufaellig orientiert)", F);
        if (kappa != 0) {                                                        // Mittel ueber die Orientierungen aus den Leistungen je Helizitaet
            real pp = 0, pm = 0; for (int a = 0; a < 3; ++a) { pp += R[a].rad_plus; pm += R[a].rad_minus; }
            const real g = 2 * (pp - pm) / (pp + pm), g0 = R[0].glum_free;
            std::printf("   g_lum gemittelt %+.5f (frei %+.5f, Verhaeltnis %.3f)", g, g0, g / g0);
        }
        std::printf("\n"); std::fflush(stdout);
    }
    return 0;
}
