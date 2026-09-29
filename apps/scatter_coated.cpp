// Streuung an beschichteten Koerpern (Kern-Schale, Mehrfachschichten): LayeredScatteringProblem.
// Kugel (Ikosaeder, Kernradius 1): Schalenflaechen als konzentrische Kugeln (Standard) oder als Parallelflaechen
// (--offset); Vergleich mit tools/mie_coated.py. Gmsh-Datei: alle Koerper erhalten dieselben Schichten.
// Schichten von innen nach aussen: --coat "d,eps_re,eps_im[,chi];..." (d in Laengeneinheiten, chi: Pasteur-Parameter);
// Kern --core "eps_re,eps_im[,chi]". --cd: beide Helizitaeten rechnen (Zeilen mit pol = +1 / -1 in der CSV, CD = Differenz).
// --inward: Schichten liegen innerhalb der angegebenen Flaeche (verdraengen Kernmaterial, z. B. Oxidation).
// Ausgabe: Q_ext bzw. sigma_ext und die Vorwaertsamplitude S(0) mit Betrag und Phase.
// --thin: Duennschicht-Naeherung erster Ordnung auf einer Flaeche (ThinLayerScatteringProblem, nur Kugel); die
//   Schichtwirkung ist dann die Differenz zu --bare (gleiches Netz), ohne neutrale Rechnung.
//   --thin-ref f: Referenzflaeche im Anteil f der Schicht von innen (0 = Innenrand, 1/2 = Mitte, 1 = Aussenrand);
//   Standard: Kernoberflaeche (nach aussen) bzw. Aussenflaeche (--inward). --bare rechnet die Kugel mit Radius 1.
//   --thin-model jump|dirac1|dirac2|dirac2fit: Sprungform 1. Ordnung, Dirac-Form 1. bzw. 2. Ordnung (Standard dirac2fit).
//   Mit --mesh und --thin: alle Koerper der Gmsh-Datei mit denselben Schichten (Referenz = Netzflaeche bzw. --thin-ref).
// --twoport: Schichten als Zweitore (S-Matrix-Formulierung, v0.19/v0.20; Schichten nach aussen, auch mehrere): E_2 auf der
//   Aussenflaeche, E_1 auf dem Kern, stabil fuer beliebiges d/h. Kugel: konzentrische Flaechen; --mesh: alle Koerper (v0.22),
//   Schichtflaechen als Parallelflaechen der Netzflaechen (Kernoberflaechen).
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
#include "cbem/problems/thin_layer_problem.hpp"
#include "cbem/problems/twoport_layer_problem.hpp"
#include "cbem/sources/fields.hpp"
using namespace cbem;
static std::vector<std::string> split(const std::string& s, char c) { std::vector<std::string> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, c)) v.push_back(t); return v; }
static cplx cval(const std::string& s) { auto c = s.find(','); return {std::stod(s.substr(0, c)), c == std::string::npos ? 0.0 : std::stod(s.substr(c + 1))}; }
int main(int argc, char** argv) {
    std::string ns = "4,6,8", path, core = "-11,1.2", coat = "0.05,2.25,0", pol = "lin", csv, nbs = "1";
    double om = 0.5, heps = 1e-4, tol = 1e-6, scale = 1.0; bool offset = false, inward = false, bare = false, neutral = false, oldnear = false, thin = false, cd = false, twoport = false; double host_chi = 0; double thinref = -1; ThinLayerModel tmodel = ThinLayerModel::Dirac2Fit; std::string tmname = "dirac2fit"; Vec3 d(0, 0, 1);
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--n") ns = nxt(); else if (o == "--mesh") path = nxt(); else if (o == "--scale") scale = std::stod(nxt());
        else if (o == "--omega") om = std::stod(nxt()); else if (o == "--core") core = nxt(); else if (o == "--coat") coat = nxt();
        else if (o == "--nbg") nbs = nxt(); else if (o == "--pol") pol = nxt(); else if (o == "--host-chi") host_chi = std::stod(nxt());
        else if (o == "--dir") { auto v = split(nxt(), ','); d = Vec3(std::stod(v[0]), std::stod(v[1]), std::stod(v[2])); d = d / norm(d); }
        else if (o == "--offset") offset = true; else if (o == "--inward") inward = true;
        else if (o == "--bare") bare = true;              // zusaetzlich dieselbe Flaeche ohne Schicht (Differenzen; Q auf denselben Radius bezogen)
        else if (o == "--thin-model") { thin = true; tmname = nxt(); tmodel = tmname == "dirac2fit" ? ThinLayerModel::Dirac2Fit : tmname == "dirac2" ? ThinLayerModel::Dirac2 : tmname == "dirac1" ? ThinLayerModel::Dirac1 : ThinLayerModel::Jump1; }
        else if (o == "--thin") thin = true; else if (o == "--cd") cd = true; else if (o == "--twoport") twoport = true; else if (o == "--thin-ref") { thin = true; thinref = std::stod(nxt()); }
        else if (o == "--neutral") neutral = true;        // zusaetzlich dieselben Netze mit Schichten aus Aussenmedium (Referenz fuer Differenzen)
        else if (o == "--oldnear") oldnear = true;        // Nahfeldregel ohne Randabstand (Vergleich der Kosten)
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt()); else if (o == "--csv") csv = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    const real nbg = std::stod(nbs); const Medium ext{nbg * nbg, 1.0, host_chi};   // --host-chi: chirales Aussenmedium (--pol circ oder --cd)
    const auto cv = split(core, ','); const Medium mcore{cval(core), 1.0, cv.size() > 2 ? std::stod(cv[2]) : 0.0};
    std::vector<Coating> coats;
    if (!coat.empty() && coat != "none")
        for (auto& c : split(coat, ';')) { auto v = split(c, ','); coats.push_back(Coating{std::stod(v[0]), Medium{cplx(std::stod(v[1]), v.size() > 2 ? std::stod(v[2]) : 0.0), 1.0, v.size() > 3 ? std::stod(v[3]) : 0.0}}); }
    real total = 0; for (auto& c : coats) total += c.thickness;
    CVec3 p = pol == "circ" ? circular_polarization(d, +1) : pol == "circ-" ? circular_polarization(d, -1) : CVec3{1.0, 0.0, 0.0};
    int cur_pol = pol == "circ" ? 1 : pol == "circ-" ? -1 : 0;
    HMatrixParams hp; hp.eps = heps; SolveOptions so; so.tol = tol;
    EntryParams ep = layered_entry_params(); if (oldnear) ep.adapt_to_boundary = false;
    bool polcol = true, chicol = true;                           // neue Spalten nur in neuen Dateien bzw. wenn vorhanden
    if (!csv.empty()) { std::ifstream in(csv); std::string h; if (std::getline(in, h)) { polcol = h.find(",pol") != std::string::npos; chicol = h.find(",core_chi") != std::string::npos; } }
    std::ofstream f; if (!csv.empty()) { f.open(csv, std::ios::app); f.seekp(0, std::ios::end); f.precision(10);
        if (f.tellp() == 0) f << "geometry,n,N,omega,core_re,core_im,coat,inward,offset,sigma_ext,Q_ext,S_re,S_im,S_abs,S_arg,iterations,near_pairs,t_near_s,t_build_s,t_solve_s" << (polcol ? ",pol" : "") << (chicol ? ",core_chi" : "") << "\n"; }
    std::printf("%s, Schichten '%s' (%s), Gesamtdicke %.4g, omega %.4g, n_bg %.3f\n", path.empty() ? "Kugel" : path.c_str(), coat.c_str(),
                inward ? "nach innen" : "nach aussen", total, om, nbg);
    std::printf("%7s %7s %11s %10s %11s %11s %9s %9s %5s %9s %8s %8s\n", "n", "N", "sigma", "Q_ext", "Re S", "Im S", "|S|", "arg S", "It.", "Nahpaare", "Nah s", "ges. s");
    auto report = [&](const std::string& label, int n, std::size_t N, const LayeredResult& r, real Aref, std::size_t np, double tn, double tb, double ts) {
        const real Q = r.sigma_ext / Aref;
        std::printf("%7s %7zu %11.6g %10.6f %11.6f %11.6f %9.6f %9.6f %5d %9zu %8.1f %8.1f\n", label.c_str(), N, r.sigma_ext, Q, r.forward.real(), r.forward.imag(),
                    std::abs(r.forward), std::arg(r.forward), r.iterations, np, tn, tb + ts);
        std::fflush(stdout);
        if (f) { f << (path.empty() ? "sphere" : path) << ',' << n << ',' << N << ',' << om << ',' << mcore.eps.real() << ',' << mcore.eps.imag() << ",\"" << (label == "bare" ? "none" : label == "neutral" ? "neutral" : label == "twoport" ? "twoport:" + coat : label == "thin" ? "thin" + (tmname != "jump" ? "-" + tmname : std::string()) + (thinref >= 0 ? "@" + std::to_string(thinref).substr(0, 4) : std::string()) + ":" + coat : coat) << "\","
                   << inward << ',' << offset << ',' << r.sigma_ext << ',' << Q << ',' << r.forward.real() << ',' << r.forward.imag() << ',' << std::abs(r.forward) << ','
                   << std::arg(r.forward) << ',' << r.iterations << ',' << np << ',' << tn << ',' << tb << ',' << ts;
                 if (polcol) f << ',' << cur_pol;
                 if (chicol) f << ',' << mcore.chi.real();
                 f << '\n'; f.flush(); }
    };
    // --cd: beide Helizitaeten auf demselben aufgebauten Problem loesen
    auto for_pols = [&](auto&& solve_one) {
        if (!cd) { solve_one(); return; }
        real sg[2]; int i = 0;
        for (int sp : {+1, -1}) { p = circular_polarization(d, sp); cur_pol = sp; sg[i++] = solve_one(); }
        std::printf("        CD = sigma+ - sigma- = %+.6e\n", sg[0] - sg[1]);
    };
    auto run_thin = [&](const std::string& label, int n, const std::vector<Coating>& cs, real Aref) {
        auto t0 = std::chrono::steady_clock::now();
        const real f = cs.empty() ? 0.0 : thinref >= 0 ? thinref : (inward ? 1.0 : 0.0);
        const real rlo = inward ? 1.0 - total : 1.0, rref = cs.empty() ? 1.0 : rlo + f * total;   // Radius der Referenzflaeche
        ThinLayerScatteringProblem P(make_icosphere(n, rref), mcore, cs, om, ext, f, hp, EntryParams{}, tmodel);
        double tb = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        for_pols([&]() { auto t1 = std::chrono::steady_clock::now(); auto r = P.solve_plane_wave(d, p, so);
            report(label, n, P.mesh().size(), r, Aref, 0, 0.0, tb, std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count()); return r.sigma_ext; });
    };
    auto run = [&](const std::string& label, int n, const LayeredGeometry& g, real Aref) {
        auto t0 = std::chrono::steady_clock::now();
        LayeredScatteringProblem P(g, om, hp, ep);
        double tb = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        for_pols([&]() { auto t1 = std::chrono::steady_clock::now(); auto r = P.solve_plane_wave(d, p, so);
            report(label, n, P.mesh().size(), r, Aref, P.near_pairs(), P.near_seconds(), tb, std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count()); return r.sigma_ext; });
    };
    std::vector<Coating> ncoats = coats; for (auto& c : ncoats) c.medium = ext;
    if (!path.empty() && twoport) {                                 // Gmsh-Koerper (alle) mit Zweitor-Schichten
        auto bodies = read_gmsh(path, scale);
        if (coats.empty() || inward) { std::printf("--twoport: Schichten nach aussen\n"); return 1; }
        auto t0 = std::chrono::steady_clock::now();
        std::vector<TwoPortBody> tb;
        for (auto& b : bodies) tb.push_back(TwoPortBody{TwoPortLayerProblem::layer_surfaces(b.mesh, coats), mcore, coats});
        TwoPortLayerProblem P(tb, om, ext, hp);
        const double tb_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        for_pols([&]() { auto t1 = std::chrono::steady_clock::now(); auto r = P.solve_plane_wave(d, p, so);
            report("twoport", 0, (coats.size() + 1) * P.outer_mesh().size(), r, 1.0, 0, 0.0, tb_s, std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count()); return r.sigma_ext; });
        return 0;
    }
    if (!path.empty() && thin) {                                    // Gmsh-Koerper mit Duennschicht-Naeherung
        auto bodies = read_gmsh(path, scale);
        const real f = thinref >= 0 ? thinref : (inward ? 1.0 : 0.0);
        const real off = inward ? -(1.0 - f) * total : f * total;
        auto solve = [&](const std::string& label, const std::vector<Coating>& cs) {
            std::vector<ThinBody> tb;
            for (auto& b : bodies) tb.push_back(ThinBody{cs.empty() || off == 0.0 ? b.mesh : offset_surface(b.mesh, off), mcore, cs, cs.empty() ? 0.0 : f});
            auto t0 = std::chrono::steady_clock::now();
            ThinLayerScatteringProblem P(tb, om, ext, hp, EntryParams{}, tmodel);
            double tb_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            for_pols([&]() { auto t1 = std::chrono::steady_clock::now(); auto r = P.solve_plane_wave(d, p, so);
                report(label, 0, P.mesh().size(), r, 1.0, 0, 0.0, tb_s, std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count()); return r.sigma_ext; });
        };
        solve("thin", coats);
        if (bare) solve("bare", {});
        return 0;
    }
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
        if (twoport) {                                               // Zweitor: Kern Radius 1, Aussenflaeche Radius 1 + d
            if (coats.empty() || inward) { std::printf("--twoport: Schichten nach aussen\n"); return 1; }
            auto t0 = std::chrono::steady_clock::now();
            std::vector<TriangleMesh> S{make_icosphere(n)}; real rr = 1.0;
            for (auto& c : coats) { rr += c.thickness; S.push_back(make_icosphere(n, rr)); }
            TwoPortLayerProblem P(S, mcore, coats, om, ext, hp);
            const double tb = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            for_pols([&]() { auto t1 = std::chrono::steady_clock::now(); auto r = P.solve_plane_wave(d, p, so);
                report("twoport", n, (coats.size() + 1) * P.outer_mesh().size(), r, pi * R * R, 0, 0.0, tb, std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count()); return r.sigma_ext; });
            if (bare) run_thin("bare", n, {}, pi * R * R);
            continue;
        }
        if (thin) {                                                  // eine Flaeche (Radius 1), Schicht erster Ordnung
            run_thin("thin", n, coats, pi * R * R);
            if (bare) run_thin("bare", n, {}, pi * R * R);
            continue;
        }
        run(std::to_string(n), n, sphere(n, coats), pi * R * R);
        if (neutral) run("neutral", n, sphere(n, ncoats), pi * R * R);
        if (bare) { LayeredGeometry g0(ext); add_body(g0, make_icosphere(n), mcore); run("bare", n, g0, pi * R * R); }
    }
    return 0;
}
