// Streuung an einem Kugel-Dimer: aeusserer Operator auf der Vereinigung aller Raender, innerer Operator
// blockdiagonal je Koerper (eigenes Medium, auch chiral); --union: ein gemeinsamer Innenoperator (gleiche
// achirale Medien, unabhaengige Kontrolle). Kugelradius 1, Mittelpunktabstand D, Einfall entlang z.
// Beispiel: scatter_multi --n 8 --dist 3 --axis x --eps1 -11,1.2 --eps2 -11,1.2 --omega 0.5 --pol circ
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/fields.hpp"
using namespace cbem;
static cplx cval(const std::string& s) { auto c = s.find(','); return {std::stod(s.substr(0, c)), c == std::string::npos ? 0.0 : std::stod(s.substr(c + 1))}; }
static std::vector<double> dlist(const std::string& s) { std::vector<double> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, ',')) v.push_back(std::stod(t)); return v; }
int main(int argc, char** argv) {
    std::string precond = "point";   // point | cluster:G (Bloecke auf Clustern mit <= G Dreiecken) | hodlr:eps[:leaf] (hierarchische Faktorisierung)
    int n = 6; std::vector<double> dists = {3.0}; std::string axis = "x", pol = "lin", csv; bool uni = false;
    double om = 0.5, heps = 1e-4, tol = 1e-6; Medium b1{cplx(-11, 1.2), 1.0, 0.0}, b2 = b1;
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--precond") precond = nxt();
        else if (o == "--n") n = std::stoi(nxt()); else if (o == "--dist") dists = dlist(nxt()); else if (o == "--axis") axis = nxt();
        else if (o == "--omega") om = std::stod(nxt()); else if (o == "--eps1") b1.eps = cval(nxt()); else if (o == "--eps2") b2.eps = cval(nxt());
        else if (o == "--chi1") b1.chi = cval(nxt()); else if (o == "--chi2") b2.chi = cval(nxt());
        else if (o == "--pol") pol = nxt(); else if (o == "--union") uni = true;
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt()); else if (o == "--csv") csv = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    std::ofstream f; if (!csv.empty()) { f.open(csv, std::ios::app); f.seekp(0, std::ios::end); if (f.tellp() == 0) f << "n,N,dist,axis,omega,eps1_re,eps1_im,chi1_re,eps2_re,eps2_im,chi2_re,pol,union,Q,Q_plus,Q_minus,iterations,MB_H,t_s\n"; }
    const TriangleMesh sphere = make_icosphere(n); const Vec3 d(0, 0, 1);
    HMatrixParams hp; hp.eps = heps; SolveOptions so; so.tol = tol;
    for (double D : dists) {
        auto t0 = std::chrono::steady_clock::now();
        Vec3 c = axis == "z" ? Vec3(0, 0, D / 2) : Vec3(D / 2, 0, 0);
        ScatteringProblem P({translated(sphere, c * -1.0), translated(sphere, c)}, {b1, b2}, om, {}, hp, {}, uni);
        if (precond.rfind("hodlr:", 0) == 0) { auto q = precond.substr(6); auto c = q.find(':'); HodlrParams hpar; hpar.eps = std::stod(q.substr(0, c)); if (c != std::string::npos) hpar.leaf = std::stoul(q.substr(c + 1)); P.use_hodlr_preconditioner(hpar); std::printf("HODLR: eps %.0e, max. Rang %zu, %.0f MB, Aufbau %.1f s\n", hpar.eps, P.hodlr()->max_rank(), P.hodlr()->bytes() / 1048576.0, P.hodlr()->seconds()); }
        else if (precond.rfind("cluster:", 0) == 0) P.use_block_preconditioner(group_by_clusters(P.mesh(), std::stoul(precond.substr(8))));
        real Q = 0, Qp = 0, Qm = 0; int its = 0;
        if (pol == "circ") {
            auto rp = P.solve_plane_wave(d, circular_polarization(d, +1), so), rm = P.solve_plane_wave(d, circular_polarization(d, -1), so);
            Qp = rp.sigma_ext / pi; Qm = rm.sigma_ext / pi; Q = 0.5 * (Qp + Qm); its = rp.iterations;
        } else { auto r = P.solve_plane_wave(d, CVec3{1.0, 0.0, 0.0}, so); Q = r.sigma_ext / pi; its = r.iterations; }
        double ts = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("D=%5.2f (%s-Achse), N=%zu: sigma_ext/pi = %.6f", D, axis.c_str(), P.mesh().size(), Q);
        if (pol == "circ") std::printf(" (+ %.6f, - %.6f, CD %.6f)", Qp, Qm, Qp - Qm);
        std::printf("  GMRES %d, H %.0f MB, %.0f s\n", its, P.hmatrix_bytes() / 1048576.0, ts); std::fflush(stdout);
        if (f) { f << n << ',' << P.mesh().size() << ',' << D << ',' << axis << ',' << om << ',' << b1.eps.real() << ',' << b1.eps.imag() << ',' << b1.chi.real() << ','
                   << b2.eps.real() << ',' << b2.eps.imag() << ',' << b2.chi.real() << ',' << pol << ',' << uni << ',' << Q << ',' << Qp << ',' << Qm << ',' << its << ','
                   << P.hmatrix_bytes() / 1048576.0 << ',' << ts << '\n'; f.flush(); }
    }
    return 0;
}
