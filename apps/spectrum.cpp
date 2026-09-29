// Spektren: Extinktion (und Zirkulardichroismus) ueber der Vakuumwellenlaenge, mit dispersiven Materialien
// (Johnson-Christy Au/Ag) in einem Hintergrundmedium, optional orientierungsgemittelt (Lebedev 6/14/26).
// Laengeneinheit der Geometrie: --unit nm (Kugel: Radius 1 Einheit; Gmsh-Datei: Koordinaten in Einheiten).
// Beispiele:
//   spectrum --sphere 12 --unit 40 --materials Au --nbg 1.33 --lambda 450:650:10 --csv results/au40_water.csv
//   spectrum --mesh examples/bornkuhn_60.msh --unit 20 --materials Au --lambda 500:900:25 --pol circ --orient 14
// Beschichtungen (alle Koerper, von innen nach aussen, Dicke in nm, Material wie --materials):
//   spectrum --sphere 12 --unit 20 --materials Ag --nbg 1.33 --coating "1.5:2.89,0" --lambda 360:460:5
//   chirale Schicht: "d:Material:chi" (Pasteur-Parameter, z. B. "1:2.25,0:0.01"); CD mit --pol circ.
//   --host-chi chi: chirales Aussenmedium (Pasteur-Parameter, z. B. chirale Loesung; nur mit --pol circ, Helizitaetswellen)
//   --twoport: Schichten als Zweitore (TwoPortLayerProblem; ein Koerper, Schichten nach aussen, beliebiges d/h).
//   --coat-inward: Schichten innerhalb der Netzflaeche (verdraengen Kernmaterial). Nur punktweise Vorkonditionierung.
//   --thin f: Duennschicht-Naeherung auf einer Flaeche je Koerper (ThinLayerScatteringProblem, auch mehrere Koerper);
//             Referenzflaeche im Anteil f der Schicht von innen (0 = Netzflaeche), per Parallelflaeche.
//   --thin-model jump|dirac1|dirac2|dirac2fit: Sprungform 1. Ordnung, Dirac-Form 1. bzw. 2. Ordnung (zusammengesetzte
//             Gradienten bzw. quadratische Anpassung; Standard dirac2fit, mit --thin 0 = Referenz auf der Metallseite).
// Ausgabe zusaetzlich: Vorwaertsamplitude S(0) (Mittel ueber Richtungen/Polarisationen) und ihre Phase arg S.
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "cbem/core/materials.hpp"
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/layered_problem.hpp"
#include "cbem/problems/thin_layer_problem.hpp"
#include "cbem/problems/twoport_layer_problem.hpp"
#include "cbem/sources/fields.hpp"
#include "cbem/sources/orientation.hpp"
using namespace cbem;
static std::vector<std::string> split(const std::string& s, char c) { std::vector<std::string> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, c)) v.push_back(t); return v; }
int main(int argc, char** argv) {
    std::string precond = "point";   // point | cluster:G (Bloecke auf Clustern mit <= G Dreiecken) | hodlr:eps[:leaf] (hierarchische Faktorisierung)
    std::string mesh, mats = "Au", chis = "0", pol = "lin", lam = "500:600:50", csv, datadir = "data/materials";
    int sph = 0, orient = 1; bool verbose = false, coat_inward = false; std::string coating; double thin = -1; bool twoport = false; double host_chi = 0; ThinLayerModel tmodel = ThinLayerModel::Dirac2Fit; std::string tmname = "dirac2fit"; double unit = 1.0, nbg = 1.0, heps = 1e-4, tol = 1e-6;
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--precond") precond = nxt();
        else if (o == "--mesh") mesh = nxt(); else if (o == "--sphere") sph = std::stoi(nxt()); else if (o == "--unit") unit = std::stod(nxt());
        else if (o == "--materials") mats = nxt(); else if (o == "--chi") chis = nxt(); else if (o == "--nbg") nbg = std::stod(nxt()); else if (o == "--host-chi") host_chi = std::stod(nxt());
        else if (o == "--lambda") lam = nxt(); else if (o == "--pol") pol = nxt(); else if (o == "--orient") orient = std::stoi(nxt());
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt());
        else if (o == "--coating") coating = nxt(); else if (o == "--coat-inward") coat_inward = true; else if (o == "--thin") thin = std::stod(nxt()); else if (o == "--twoport") twoport = true;
        else if (o == "--thin-model") { tmname = nxt(); tmodel = tmname == "jump" ? ThinLayerModel::Jump1 : tmname == "dirac1" ? ThinLayerModel::Dirac1 : tmname == "dirac2fit" ? ThinLayerModel::Dirac2Fit : ThinLayerModel::Dirac2; }
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
    std::vector<std::pair<real, std::shared_ptr<Material>>> coat_mat;           // Dicke in nm, Material
    std::vector<real> coat_chi;                                                   // Pasteur-Parameter je Schicht (optional ":chi")
    for (auto& c : split(coating, ';')) {
        if (c.empty()) continue; auto q = c.find(':'); auto q2 = c.find(':', q + 1);
        coat_mat.push_back({std::stod(c.substr(0, q)), make_material(c.substr(q + 1, q2 == std::string::npos ? std::string::npos : q2 - q - 1), datadir)});
        coat_chi.push_back(q2 == std::string::npos ? 0.0 : std::stod(c.substr(q2 + 1)));
    }
    if (!coat_mat.empty() && precond != "point") { std::printf("Beschichtung: nur --precond point\n"); return 1; }
    if (thin >= 0 && coat_mat.empty()) { std::printf("--thin: nur mit --coating\n"); return 1; }
    if (twoport && (coat_mat.empty() || coat_inward)) { std::printf("--twoport: nur mit --coating nach aussen\n"); return 1; }
    std::vector<TriangleMesh> tp_surf;                         // Zweitor: Kern und Schichtflaechen (Parallelflaechen, einmal)
    if (twoport) { std::vector<Coating> cs; for (auto& c : coat_mat) cs.push_back(Coating{c.first / unit, Medium{}}); tp_surf = TwoPortLayerProblem::layer_surfaces(parts[0], cs);
                   if (parts.size() != 1) std::printf("Hinweis: --twoport rechnet nur den ersten Koerper\n"); }
    std::vector<TriangleMesh> thin_ref;                        // Referenzflaechen der Duennschicht-Naeherung (je Koerper)
    if (thin >= 0) { real t = 0; for (auto& c : coat_mat) t += c.first / unit; const real off = coat_inward ? -(1.0 - thin) * t : thin * t;
                     for (auto& pm : parts) thin_ref.push_back(off == 0.0 ? pm : offset_surface(pm, off)); }
    std::vector<double> lams; { auto r = split(lam, ':'); if (r.size() == 3) for (double l = std::stod(r[0]); l <= std::stod(r[1]) + 1e-9; l += std::stod(r[2])) lams.push_back(l); else for (auto& t : split(lam, ',')) lams.push_back(std::stod(t)); }
    const auto dirs = lebedev(orient);
    std::ofstream f; if (!csv.empty()) { f.open(csv, std::ios::app); f.precision(10); f.seekp(0, std::ios::end); if (f.tellp() == 0) f << "lambda_nm,unit_nm,nbg,N,orient,pol,sigma_nm2,sigma_plus_nm2,sigma_minus_nm2,CD_nm2,iterations,t_s,coating,S_re,S_im\n"; }
    std::printf("%s: %zu Koerper, %zu Dreiecke gesamt; Einheit %.3g nm, n_Hintergrund %.3f, %zu Richtung(en)\n",
                sph ? "Kugel" : mesh.c_str(), parts.size(), [&] { std::size_t n = 0; for (auto& p : parts) n += p.size(); return n; }(), unit, nbg, dirs.size());
    if (!coat_mat.empty()) std::printf("Beschichtung '%s' (%s)\n", coating.c_str(), coat_inward ? "nach innen" : "nach aussen");
    std::printf("%9s %13s %13s %13s %13s %9s %6s %7s\n", "lambda nm", "sigma nm^2", "sigma+ nm^2", "sigma- nm^2", "CD nm^2", "arg S", "It.", "Zeit s");
    for (double L : lams) {
        auto t0 = std::chrono::steady_clock::now();
        const real om = 2 * pi * unit / L;                     // k0 in 1/Einheit
        std::vector<Medium> med;
        for (std::size_t b = 0; b < parts.size(); ++b) med.push_back(Medium{mat[b]->eps(L), 1.0, chi[b]});
        HMatrixParams hp; hp.eps = heps; SolveOptions so; so.tol = tol;
        const Medium bg{nbg * nbg, 1.0, host_chi};                      // --host-chi: chirales Aussenmedium (nur --pol circ)
        std::unique_ptr<ScatteringProblem> PP; std::unique_ptr<LayeredScatteringProblem> PL; std::unique_ptr<ThinLayerScatteringProblem> PT;
        std::unique_ptr<TwoPortLayerProblem> P2;
        if (coat_mat.empty()) PP = std::make_unique<ScatteringProblem>(parts, med, om, bg, hp);
        else if (twoport) {
            std::vector<Coating> cs; for (std::size_t i = 0; i < coat_mat.size(); ++i) cs.push_back(Coating{coat_mat[i].first / unit, Medium{coat_mat[i].second->eps(L), 1.0, coat_chi[i]}});
            P2 = std::make_unique<TwoPortLayerProblem>(tp_surf, med[0], cs, om, bg, hp);
        }
        else if (thin >= 0) {
            std::vector<Coating> cs; for (std::size_t i = 0; i < coat_mat.size(); ++i) cs.push_back(Coating{coat_mat[i].first / unit, Medium{coat_mat[i].second->eps(L), 1.0, coat_chi[i]}});
            std::vector<ThinBody> tb;
            for (std::size_t b = 0; b < parts.size(); ++b) tb.push_back(ThinBody{thin_ref[b], med[b], cs, thin});
            PT = std::make_unique<ThinLayerScatteringProblem>(tb, om, bg, hp, EntryParams{}, tmodel);
        }
        else {
            LayeredGeometry g(bg); std::vector<Coating> cs;
            for (std::size_t i = 0; i < coat_mat.size(); ++i) cs.push_back(Coating{coat_mat[i].first / unit, Medium{coat_mat[i].second->eps(L), 1.0, coat_chi[i]}});
            for (std::size_t b = 0; b < parts.size(); ++b) add_coated_body(g, parts[b], med[b], cs, !coat_inward);
            PL = std::make_unique<LayeredScatteringProblem>(g, om, hp);
        }
        const std::size_t Ntri = PP ? PP->mesh().size() : PT ? PT->mesh().size() : P2 ? (P2->layers() + 1) * P2->outer_mesh().size() : PL->mesh().size();
        struct Sol { real sigma_ext; cplx forward; int iterations; };
        auto solve = [&](const Vec3& dd, const CVec3& pp) -> Sol {
            if (PP) { auto r = PP->solve_plane_wave(dd, pp, so); return {r.sigma_ext, r.forward, r.iterations}; }
            if (PT) { auto r = PT->solve_plane_wave(dd, pp, so); return {r.sigma_ext, r.forward, r.iterations}; }
            if (P2) { auto r = P2->solve_plane_wave(dd, pp, so); return {r.sigma_ext, r.forward, r.iterations}; }
            auto r = PL->solve_plane_wave(dd, pp, so); return {r.sigma_ext, r.forward, r.iterations}; };
        if (PP) {
            ScatteringProblem& P = *PP;
            if (precond.rfind("hodlr:", 0) == 0) { auto q = precond.substr(6); auto c = q.find(':'); HodlrParams hpar; hpar.eps = std::stod(q.substr(0, c)); if (c != std::string::npos) hpar.leaf = std::stoul(q.substr(c + 1)); P.use_hodlr_preconditioner(hpar); std::printf("HODLR: eps %.0e, max. Rang %zu, %.0f MB, Aufbau %.1f s\n", hpar.eps, P.hodlr()->max_rank(), P.hodlr()->bytes() / 1048576.0, P.hodlr()->seconds()); }
            else if (precond.rfind("cluster:", 0) == 0) P.use_block_preconditioner(group_by_clusters(P.mesh(), std::stoul(precond.substr(8))));
        }
        cplx Sf = 0;
        real s = 0, sp = 0, sm = 0; int its = 0;
        for (const auto& dw : dirs) {
            if (pol == "circ") {
                auto a = solve(dw.d, circular_polarization(dw.d, +1)), b = solve(dw.d, circular_polarization(dw.d, -1));
                sp += dw.w * a.sigma_ext; sm += dw.w * b.sigma_ext; Sf += dw.w * 0.5 * (a.forward + b.forward); its = std::max(its, std::max(a.iterations, b.iterations));
                if (verbose) std::printf("      d = (%+.3f, %+.3f, %+.3f): sigma+ %.5g, sigma- %.5g, CD %+.4g nm^2\n", dw.d.x, dw.d.y, dw.d.z,
                                         a.sigma_ext * unit * unit, b.sigma_ext * unit * unit, (a.sigma_ext - b.sigma_ext) * unit * unit);
            } else {
                // linear: Mittel ueber zwei orthogonale Polarisationen = unpolarisiert (fuer orient > 1), sonst p = u
                CVec3 u = circular_polarization(dw.d, +1), v = circular_polarization(dw.d, -1);
                CVec3 pu{0.5 * (u[0] + v[0]), 0.5 * (u[1] + v[1]), 0.5 * (u[2] + v[2])};
                auto a = solve(dw.d, pu); real sa = a.sigma_ext; cplx Sa = a.forward; its = std::max(its, a.iterations);
                if (dirs.size() > 1) {
                    CVec3 pv{(u[0] - v[0]) / cplx(0, 2), (u[1] - v[1]) / cplx(0, 2), (u[2] - v[2]) / cplx(0, 2)};
                    auto bv = solve(dw.d, pv); sa = 0.5 * (sa + bv.sigma_ext); Sa = 0.5 * (Sa + bv.forward);
                }
                Sf += dw.w * Sa;
                s += dw.w * sa;
            }
        }
        if (pol == "circ") s = 0.5 * (sp + sm);
        const real u2 = unit * unit;
        double ts = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("%9.1f %13.5g %13.5g %13.5g %13.5g %9.5f %6d %7.1f\n", L, s * u2, sp * u2, sm * u2, (sp - sm) * u2, std::arg(Sf), its, ts); std::fflush(stdout);
        if (f) { f << L << ',' << unit << ',' << nbg << ',' << Ntri << ',' << dirs.size() << ',' << pol << ',' << s * u2 << ',' << sp * u2 << ',' << sm * u2 << ','
                   << (sp - sm) * u2 << ',' << its << ',' << ts << ",\"" << coating << (coat_inward ? " (innen)" : "") << (thin >= 0 ? " thin-" + tmname + "@" + std::to_string(thin).substr(0, 4) : std::string()) << (twoport ? " twoport" : "") << "\"," << Sf.real() << ',' << Sf.imag() << '\n'; f.flush(); }
    }
    return 0;
}
