// Streuung an beschichteten Koerpern (Kern-Schale, Mehrfachschichten): LayeredScatteringProblem.
// Kugel (Ikosaeder, Kernradius 1): Schalenflaechen als konzentrische Kugeln (Standard) oder als Parallelflaechen
// (--offset); Vergleich mit tools/mie_coated.py. Gmsh-Datei: alle Koerper erhalten dieselben Schichten.
// Schichten von innen nach aussen: --coat "d,eps_re,eps_im;d,eps_re,eps_im" (d in Laengeneinheiten).
// --inward: Schichten liegen innerhalb der angegebenen Flaeche (verdraengen Kernmaterial, z. B. Oxidation).
// Ausgabe: Q_ext bzw. sigma_ext und die Vorwaertsamplitude S(0) mit Betrag und Phase.
// --bare: dieselbe Flaeche ohne Schicht; --neutral: dieselben Netze, Schichten aus dem Aussenmedium. Die Wirkung einer
// duennen Schicht ist als Differenz zur neutralen Rechnung genauer als zur Rechnung ohne Schicht (docs/results_coated.md).
// Beispiele:
//   scatter_coated --n 4,8,12 --omega 0.5 --core -11,1.2 --coat 0.05,2.25,0
//   scatter_coated --mesh roundcube.msh --omega 0.8 --core -8,0.5 --coat 0.04,3.0,0 --inward --dir 0,0,1
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/layered_problem.hpp"
#include "cbem/sources/fields.hpp"
using namespace cbem;
static std::vector<std::string> split(const std::string& s, char c) { std::vector<std::string> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, c)) v.push_back(t); return v; }
static cplx cval(const std::string& s) { auto c = s.find(','); return {std::stod(s.substr(0, c)), c == std::string::npos ? 0.0 : std::stod(s.substr(c + 1))}; }
int main(int argc, char** argv) {
    std::string ns = "4,6,8", path, core = "-11,1.2", coat = "0.05,2.25,0", pol = "lin", csv, nbs = "1";
    double om = 0.5, heps = 1e-4, tol = 1e-6, scale = 1.0; bool offset = false, inward = false, bare = false, neutral = false, oldnear = false; Vec3 d(0, 0, 1);
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--n") ns = nxt(); else if (o == "--mesh") path = nxt(); else if (o == "--scale") scale = std::stod(nxt());
        else if (o == "--omega") om = std::stod(nxt()); else if (o == "--core") core = nxt(); else if (o == "--coat") coat = nxt();
        else if (o == "--nbg") nbs = nxt(); else if (o == "--pol") pol = nxt();
        else if (o == "--dir") { auto v = split(nxt(), ','); d = Vec3(std::stod(v[0]), std::stod(v[1]), std::stod(v[2])); d = d / norm(d); }
        else if (o == "--offset") offset = true; else if (o == "--inward") inward = true;
        else if (o == "--bare") bare = true;              // zusaetzlich dieselbe Flaeche ohne Schicht (Differenzen; Q auf denselben Radius bezogen)
        else if (o == "--neutral") neutral = true;        // zusaetzlich dieselben Netze mit Schichten aus Aussenmedium (Referenz fuer Differenzen)
        else if (o == "--oldnear") oldnear = true;        // Nahfeldregel ohne Randabstand (Vergleich der Kosten)
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt()); else if (o == "--csv") csv = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    const real nbg = std::stod(nbs); const Medium ext{nbg * nbg, 1.0, 0.0}, mcore{cval(core), 1.0, 0.0};
    std::vector<Coating> coats;
    if (!coat.empty() && coat != "none")
        for (auto& c : split(coat, ';')) { auto v = split(c, ','); coats.push_back(Coating{std::stod(v[0]), Medium{cplx(std::stod(v[1]), v.size() > 2 ? std::stod(v[2]) : 0.0), 1.0, 0.0}}); }
    real total = 0; for (auto& c : coats) total += c.thickness;
    const CVec3 p = pol == "circ" ? circular_polarization(d, +1) : pol == "circ-" ? circular_polarization(d, -1) : CVec3{1.0, 0.0, 0.0};
    HMatrixParams hp; hp.eps = heps; SolveOptions so; so.tol = tol;
    EntryParams ep = layered_entry_params(); if (oldnear) ep.adapt_to_boundary = false;
    std::ofstream f; if (!csv.empty()) { f.open(csv, std::ios::app); f.seekp(0, std::ios::end);
        if (f.tellp() == 0) f << "geometry,n,N,omega,core_re,core_im,coat,inward,offset,sigma_ext,Q_ext,S_re,S_im,S_abs,S_arg,iterations,near_pairs,t_near_s,t_build_s,t_solve_s\n"; }
    std::printf("%s, Schichten '%s' (%s), Gesamtdicke %.4g, omega %.4g, n_bg %.3f\n", path.empty() ? "Kugel" : path.c_str(), coat.c_str(),
                inward ? "nach innen" : "nach aussen", total, om, nbg);
    std::printf("%7s %7s %11s %10s %11s %11s %9s %9s %5s %9s %8s %8s\n", "n", "N", "sigma", "Q_ext", "Re S", "Im S", "|S|", "arg S", "It.", "Nahpaare", "Nah s", "ges. s");
    auto run = [&](const std::string& label, int n, const LayeredGeometry& g, real Aref) {
        auto t0 = std::chrono::steady_clock::now();
        LayeredScatteringProblem P(g, om, hp, ep);
        double tb = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); t0 = std::chrono::steady_clock::now();
        auto r = P.solve_plane_wave(d, p, so);
        double ts = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const real Q = r.sigma_ext / Aref;
        std::printf("%7s %7zu %11.6g %10.6f %11.6f %11.6f %9.6f %9.6f %5d %9zu %8.1f %8.1f\n", label.c_str(), P.mesh().size(), r.sigma_ext, Q, r.forward.real(), r.forward.imag(),
                    std::abs(r.forward), std::arg(r.forward), r.iterations, P.near_pairs(), P.near_seconds(), tb + ts);
        std::fflush(stdout);
        if (f) { f << (path.empty() ? "sphere" : path) << ',' << n << ',' << P.mesh().size() << ',' << om << ',' << mcore.eps.real() << ',' << mcore.eps.imag() << ",\"" << (label == "bare" ? "none" : label == "neutral" ? "neutral" : coat) << "\","
                   << inward << ',' << offset << ',' << r.sigma_ext << ',' << Q << ',' << r.forward.real() << ',' << r.forward.imag() << ',' << std::abs(r.forward) << ','
                   << std::arg(r.forward) << ',' << r.iterations << ',' << P.near_pairs() << ',' << P.near_seconds() << ',' << tb << ',' << ts << '\n'; f.flush(); }
    };
    std::vector<Coating> ncoats = coats; for (auto& c : ncoats) c.medium = ext;
    if (!path.empty()) {
        auto bodies = read_gmsh(path, scale);
        LayeredGeometry g(ext), g0(ext), gn(ext);
        for (auto& b : bodies) { add_coated_body(g, b.mesh, mcore, coats, !inward); add_body(g0, b.mesh, mcore); add_coated_body(gn, b.mesh, mcore, ncoats, !inward); }
        run("mesh", 0, g, 1.0);
        if (bare) run("bare", 0, g0, 1.0);
        if (neutral) run("neutral", 0, gn, 1.0);
        return 0;
    }
    const real R = inward ? 1.0 : 1.0 + total;                    // Aussenradius; Q_ext bezogen auf pi R^2
    auto sphere = [&](int n, const std::vector<Coating>& cs) {
        LayeredGeometry g(ext);
        if (offset) add_coated_body(g, make_icosphere(n), mcore, cs, !inward);
        else {
            std::vector<TriangleMesh> surf; std::vector<Medium> med; real r = R;
            surf.push_back(make_icosphere(n, r));
            for (std::size_t l = cs.size(); l-- > 0;) { med.push_back(cs[l].medium); r -= cs[l].thickness; surf.push_back(make_icosphere(n, r)); }
            med.push_back(mcore);
            add_layered_body(g, surf, med);
        }
        return g;
    };
    for (auto& s : split(ns, ',')) {
        const int n = std::stoi(s);
        run(std::to_string(n), n, sphere(n, coats), pi * R * R);
        if (neutral) run("neutral", n, sphere(n, ncoats), pi * R * R);
        if (bare) { LayeredGeometry g0(ext); add_body(g0, make_icosphere(n), mcore); run("bare", n, g0, pi * R * R); }
    }
    return 0;
}
