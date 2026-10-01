// Optische Pinzette (v0.37): Kraft eines stark fokussierten Laserstrahls (Richards-Wolf, exaktes Winkelspektrum) oder eines
// Gaussstrahls auf eine Kugel. Das Teilchen bleibt fest, der Brennpunkt wandert: Der Operator wird einmal aufgebaut, jede
// Brennpunktlage ist eine neue rechte Seite. Kraft ueber den Spannungstensor auf einer Kugel um das Teilchen, ausgegeben als
// Effizienz Q = F c / (n P) (Strahlleistung P = 1). Positionen sind die Lage des Teilchens relativ zum Brennpunkt (in um).
// Beschichtete Kugeln (v0.41): --shell "t:n" (Schichtdicke t in um, Brechungsindex n; Zweitor-Formulierung); --radius ist dann
// der Kernradius. Kraft auf einer Kugel um die aeussere Flaeche.
// Chirale Medien (v0.38): --host-chi (chirales Aussenmedium, Strahl nach Helizitaeten zerlegt), --chi-particle (chirales Teilchen);
// --pol circ+ / circ- (zirkular in der Pupille).
// Beispiel: tweezers --radius 0.25 --n-particle 1.59 --nbg 1.33 --lambda 1.064 --NA 1.2 --f0 1 --axial "-1:2:13" --mesh 8
//           tweezers ... --lateral "0:0.6:7" --at-z 0.3
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/problems/twoport_layer_problem.hpp"
#include "cbem/sources/optical_force.hpp"
using namespace cbem;
static std::vector<std::string> split(const std::string& s, char c) { std::vector<std::string> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, c)) v.push_back(t); return v; }
int main(int argc, char** argv) {
    double R = 0.25, np = 1.59, nk = 0.0, nbg = 1.33, lam = 1.064, NA = 1.2, f0 = 1.0, w0 = -1, atz = 0, heps = 1e-6, tol = 1e-8, hchi = 0, pchi = 0;
    int nmesh = 8, nt = 40; std::string pol = "x", axial, lateral, csv, shell;
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--radius") R = std::stod(nxt()); else if (o == "--n-particle") np = std::stod(nxt()); else if (o == "--k-particle") nk = std::stod(nxt());
        else if (o == "--nbg") nbg = std::stod(nxt()); else if (o == "--lambda") lam = std::stod(nxt()); else if (o == "--NA") NA = std::stod(nxt());
        else if (o == "--f0") f0 = std::stod(nxt()); else if (o == "--gaussian") w0 = std::stod(nxt()); else if (o == "--pol") pol = nxt();
        else if (o == "--axial") axial = nxt(); else if (o == "--lateral") lateral = nxt(); else if (o == "--at-z") atz = std::stod(nxt());
        else if (o == "--host-chi") hchi = std::stod(nxt()); else if (o == "--chi-particle") pchi = std::stod(nxt());
        else if (o == "--shell") shell = nxt();
        else if (o == "--mesh") nmesh = std::stoi(nxt()); else if (o == "--ntheta") nt = std::stoi(nxt());
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt()); else if (o == "--csv") csv = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    // Einheiten: Laengen in Teilchenradien
    const real om = 2 * pi * R / lam;
    const Medium host{nbg * nbg, 1.0, hchi}, part{cplx(np, nk) * cplx(np, nk), 1.0, pchi};
    const real s2 = 1.0 / std::sqrt(2.0);
    const CVec3 pp = (pol == "circ" || pol == "circ+") ? CVec3{s2, cplx(0, s2), 0.0} : pol == "circ-" ? CVec3{s2, cplx(0, -s2), 0.0} : pol == "y" ? CVec3{0.0, 1.0, 0.0} : CVec3{1.0, 0.0, 0.0};
    HMatrixParams hp; hp.eps = heps; SolveOptions so; so.tol = tol;
    auto t0 = std::chrono::steady_clock::now();
    // homogen: ScatteringProblem; beschichtet: Zweitor (Kern Radius 1, Schicht t/R)
    std::unique_ptr<ScatteringProblem> PS; std::unique_ptr<TwoPortLayerProblem> P2; TriangleMesh outer; real Rout = 1.0;
    if (shell.empty()) { PS = std::make_unique<ScatteringProblem>(std::vector<TriangleMesh>{make_icosphere(nmesh)}, std::vector<Medium>{part}, om, host, hp); outer = PS->mesh(); }
    else {
        const auto v = split(shell, ':'); const real t = std::stod(v[0]) / R, ns = std::stod(v[1]); Rout = 1 + t;
        const Medium sh{ns * ns, 1.0, 0.0}; const TriangleMesh core = make_icosphere(nmesh);
        const auto surf = TwoPortLayerProblem::layer_surfaces(core, {Coating{t, sh}});
        P2 = std::make_unique<TwoPortLayerProblem>(surf, part, std::vector<Coating>{Coating{t, sh}}, om, host, hp); outer = P2->outer_mesh();
        std::printf("Schicht %.3f um mit n = %.3f, aeusserer Radius %.3f um\n", t * R, ns, Rout * R);
    }
    struct { std::size_t n; } meshsz{outer.size()};
    if (hchi != 0 || pchi != 0) std::printf("chiral: Aussenmedium chi = %g, Teilchen chi = %g, Polarisation %s\n", hchi, pchi, pol.c_str());
    std::printf("Kugel R = %.3f um, n = %.3f%+.3fi in n = %.3f, lambda = %.3f um, %s, %zu Dreiecke (k R = %.2f); Aufbau %.1f s\n", R, np, nk, nbg, lam,
                w0 > 0 ? ("Gaussstrahl w0 = " + std::to_string(w0) + " um").c_str() : ("NA " + std::to_string(NA).substr(0, 4) + ", f0 " + std::to_string(f0).substr(0, 4)).c_str(),
                meshsz.n, std::real(host.k(om)) * Rout, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    std::vector<Vec3> pos;                                                       // Teilchen relativ zum Brennpunkt (um)
    auto rng = [&](const std::string& s) { auto v = split(s, ':'); std::vector<real> out; const int n = std::stoi(v[2]);
        for (int i = 0; i < n; ++i) out.push_back(std::stod(v[0]) + (std::stod(v[1]) - std::stod(v[0])) * i / std::max(1, n - 1)); return out; };
    if (!axial.empty()) for (real z : rng(axial)) pos.push_back(Vec3(0, 0, z));
    if (!lateral.empty()) for (real x : rng(lateral)) pos.push_back(Vec3(x, 0, atz));
    std::ofstream f; if (!csv.empty()) { f.open(csv); f << "x_um,y_um,z_um,Qx,Qy,Qz,iterations\n"; }
    std::printf("      x um      z um          Qx          Qy          Qz   It.   s\n");
    for (const Vec3& x : pos) {
        t0 = std::chrono::steady_clock::now();
        const Vec3 focus = x * (-1.0 / R);                                       // Brennpunkt relativ zum Teilchen (Einheiten R)
        auto beam = std::make_shared<BeamField>(w0 > 0 ? BeamField::gaussian(host, om, focus, w0 / R, pp, nt, 2 * nt)
                                                       : BeamField::focused(host, om, focus, NA, f0, pp, nt, 2 * nt));
        const auto b = beam->project(outer, host);
        std::vector<cplx> h; int its;
        if (PS) { const auto r = PS->solve_rhs(b, so); h = r.h; its = r.iterations; }
        else { const auto r = P2->solve_rhs(b, so); h.assign(r.h.begin(), r.h.begin() + b.size()); its = r.iterations; }
        const Vec3 F = force_on_sphere(make_near_field_eval(outer, h, b, host, om, beam), host, Vec3(0, 0, 0), 1.3 * Rout);
        const Vec3 Q = F / std::real(std::sqrt(host.eps * host.mu));             // Leistung 1, c = 1
        std::printf("  %8.3f  %8.3f  %+.4e %+.4e %+.4e  %4d  %.1f\n", x.x, x.z, Q.x, Q.y, Q.z, its, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        std::fflush(stdout);
        if (f) f << x.x << ',' << x.y << ',' << x.z << ',' << Q.x << ',' << Q.y << ',' << Q.z << ',' << its << std::endl;
    }
    return 0;
}
