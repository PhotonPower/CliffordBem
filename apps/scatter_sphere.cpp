// Streuung einer ebenen Welle an der (Ikosaeder-)Kugel mit T_1, H-Matrizen und GMRES.
// Beispiel: scatter_sphere --n 8,16,24 --omega 1.0 --eps1 2.25,0 --csv results/scatter_glass.csv
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/solvers/gmres.hpp"
#include "cbem/sources/fields.hpp"
using namespace cbem;
static std::vector<int> parse_list(const std::string& s) { std::vector<int> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, ',')) v.push_back(std::stoi(t)); return v; }
int main(int argc, char** argv) {
    std::vector<int> ns = {4, 6, 8}; double om = 1.0, e1r = 2.25, e1i = 0.0, heps = 1e-4, tol = 1e-6; int restart = 200; std::string csv;
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--n") ns = parse_list(nxt()); else if (o == "--omega") om = std::stod(nxt());
        else if (o == "--eps1") { std::string s = nxt(); auto c = s.find(','); e1r = std::stod(s.substr(0, c)); e1i = c == std::string::npos ? 0 : std::stod(s.substr(c + 1)); }
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt());
        else if (o == "--restart") restart = std::stoi(nxt()); else if (o == "--csv") csv = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    const cplx eps1(e1r, e1i); Medium in{eps1, 1.0}, out{1.0, 1.0};
    const cplx k1 = om * std::sqrt(eps1), k2 = om;
    std::ofstream f; if (!csv.empty()) { f.open(csv, std::ios::app); f.seekp(0, std::ios::end); if (f.tellp() == 0) f << "n,N,unknowns,omega,eps1_re,eps1_im,heps,Qext,iterations,residual,MB_H,t_build_s,t_solve_s\n"; }
    std::printf("%5s %8s %9s %12s %6s %10s %8s %9s %9s\n", "n", "N", "Unbek.", "Q_ext", "It.", "Residuum", "MB H", "Aufbau s", "Loesen s");
    for (int n : ns) {
        TriangleMesh m = make_icosphere(n);
        auto t0 = std::chrono::steady_clock::now();
        KernelEntries Ein(m, k1), Eout(m, k2); HMatrixParams p; p.eps = heps;
        KernelHMatrix H1(Ein, p), H2(Eout, p); CauchyOperator E1(m, H1), E2(m, H2);
        TransmissionOperator T(m, E1, E2, in, out);
        double tb = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        Vec3 d(0, 0, 1), pol(1, 0, 0);
        auto b = project_plane_wave(m, k2, 1.0, d, pol);
        LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.apply(x, y); };
        LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.precondition(x, y); };
        t0 = std::chrono::steady_clock::now();
        std::vector<cplx> h; GmresResult r = gmres(A, b, h, &M, tol, restart, 3000);
        double ts = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::vector<cplx> hs(h.size()); for (std::size_t i = 0; i < h.size(); ++i) hs[i] = h[i] - b[i];
        real Q = extinction_cross_section(m, hs, k2, 1.0, d, pol) / pi;     // Kugelradius 1
        double MB = (H1.stats().bytes() + H2.stats().bytes()) / 1048576.0;
        std::printf("%5d %8zu %9zu %12.6f %6d %10.1e %8.0f %9.1f %9.1f\n", n, m.size(), 8 * m.size(), Q, r.iterations, r.rel_residual, MB, tb, ts);
        std::fflush(stdout);
        if (f) { f << n << ',' << m.size() << ',' << 8 * m.size() << ',' << om << ',' << e1r << ',' << e1i << ',' << heps << ',' << Q << ',' << r.iterations << ','
                   << r.rel_residual << ',' << MB << ',' << tb << ',' << ts << '\n'; f.flush(); }
    }
    return 0;
}
