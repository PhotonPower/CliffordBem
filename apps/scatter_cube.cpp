// Streuung am Wuerfel [-1,1]^3 mit T_1, H-Matrizen (ACA mit exakten Eintraegen, ohne h-Kriterium)
// und Vergleich punktweiser gegen Kanten-/Eck-Blockvorkonditionierung.
// Beispiel: scatter_cube --mesh graded --L 3,5,7 --omega 0.5 --eps1 -11,1.2 --R 0.25 --csv results/cube_gold.csv
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/solvers/block_preconditioner.hpp"
#include "cbem/solvers/gmres.hpp"
#include "cbem/sources/fields.hpp"
using namespace cbem;
static std::vector<int> ilist(const std::string& s) { std::vector<int> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, ',')) v.push_back(std::stoi(t)); return v; }
static std::vector<double> dlist(const std::string& s) { std::vector<double> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, ',')) v.push_back(std::stod(t)); return v; }
int main(int argc, char** argv) {
    std::string meshkind = "graded", csv; std::vector<int> Ls = {3, 5}; std::vector<double> Rs = {0.25};
    double om = 0.5, e1r = -11, e1i = 1.2, heps = 1e-4, tol = 1e-6; bool point = true;
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--mesh") meshkind = nxt(); else if (o == "--L") Ls = ilist(nxt()); else if (o == "--R") Rs = dlist(nxt());
        else if (o == "--omega") om = std::stod(nxt());
        else if (o == "--eps1") { std::string s = nxt(); auto c = s.find(','); e1r = std::stod(s.substr(0, c)); e1i = c == std::string::npos ? 0 : std::stod(s.substr(c + 1)); }
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt());
        else if (o == "--no-point") point = false; else if (o == "--csv") csv = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    const cplx eps1(e1r, e1i); Medium in{eps1, 1.0}, out{1.0, 1.0}; const cplx k1 = om * std::sqrt(eps1), k2 = om;
    std::ofstream f; if (!csv.empty()) { f.open(csv, std::ios::app); f.seekp(0, std::ios::end); if (f.tellp() == 0) f << "mesh,L,N,omega,eps1_re,eps1_im,precond,R,max_block,block_fraction,t_setup_s,iterations,residual,sigma_ext,MB_H,t_H_s\n"; }
    const Vec3 d(0, 0, 1), pol(1, 0, 0);
    for (int L : Ls) {
        TriangleMesh m = meshkind == "uniform" ? make_cube_uniform(L) : make_cube_graded(L);
        auto t0 = std::chrono::steady_clock::now();
        KernelEntries Ein(m, k1), Eout(m, k2);
        HMatrixParams p; p.eps = heps; p.sep_factor = 0.0; p.exact_in_lowrank = true;
        KernelHMatrix H1(Ein, p), H2(Eout, p); CauchyOperator E1(m, H1), E2(m, H2);
        TransmissionOperator T(m, E1, E2, in, out);
        double tH = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        double MB = (H1.stats().bytes() + H2.stats().bytes()) / 1048576.0;
        auto b = project_plane_wave(m, k2, 1.0, d, pol);
        LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.apply(x, y); };
        std::printf("%s L=%d: N=%zu (%zu Unbekannte), H %.0f MB (%.0f s)\n", meshkind.c_str(), L, m.size(), 8 * m.size(), MB, tH);
        auto run = [&](const LinOp& M, const std::string& name, double R, std::size_t mb, double fr, double ts) {
            std::vector<cplx> h; GmresResult r = gmres(A, b, h, &M, tol, 300, 3000);
            std::vector<cplx> hs(h.size()); for (std::size_t i = 0; i < h.size(); ++i) hs[i] = h[i] - b[i];
            real sig = extinction_cross_section(m, hs, k2, 1.0, d, pol);
            std::printf("   %-9s R=%.3f: GMRES %4d It. (Residuum %.1e), sigma_ext = %.6f, groesster Block %zu, Blockanteil %.2f, Aufbau %.1f s\n",
                        name.c_str(), R, r.iterations, r.rel_residual, sig, mb, fr, ts);
            std::fflush(stdout);
            if (f) { f << meshkind << ',' << L << ',' << m.size() << ',' << om << ',' << e1r << ',' << e1i << ',' << name << ',' << R << ',' << mb << ',' << fr << ','
                       << ts << ',' << r.iterations << ',' << r.rel_residual << ',' << sig << ',' << MB << ',' << tH << '\n'; f.flush(); }
        };
        if (point) { LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.precondition(x, y); }; run(M, "punktweise", 0.0, 8, 0.0, 0.0); }
        for (double R : Rs) {
            auto groups = group_by_features(m, FeatureSet::cube(1.0), R);
            BlockPreconditioner P(m, Ein, Eout, T, groups);
            LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { P.apply(x, y); };
            run(M, "Bloecke", R, P.max_block(), P.fraction(), P.seconds());
        }
    }
    return 0;
}
