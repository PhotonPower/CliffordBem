// Zirkulardichroismus einer chiralen Kugel: sigma_ext fuer beide zirkularen Polarisationen, CD = sigma_+ - sigma_-.
// Beispiel: scatter_chiral --n 8,12,16 --omega 1.0 --eps1 2.25,0 --chi 0.2,0 --csv results/chiral_glass.csv
#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "cbem/operators/chiral_cauchy_operator.hpp"
#include "cbem/operators/transmission_operator.hpp"
#include "cbem/solvers/gmres.hpp"
#include "cbem/sources/fields.hpp"
using namespace cbem;
static std::vector<int> ilist(const std::string& s) { std::vector<int> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, ',')) v.push_back(std::stoi(t)); return v; }
static cplx cval(const std::string& s) { auto c = s.find(','); return {std::stod(s.substr(0, c)), c == std::string::npos ? 0.0 : std::stod(s.substr(c + 1))}; }
int main(int argc, char** argv) {
    std::vector<int> ns = {4, 6, 8}; double om = 1.0, heps = 1e-4, tol = 1e-6; cplx eps1 = 2.25, chi1 = 0.2; std::string csv;
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--n") ns = ilist(nxt()); else if (o == "--omega") om = std::stod(nxt());
        else if (o == "--eps1") eps1 = cval(nxt()); else if (o == "--chi") chi1 = cval(nxt());
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt()); else if (o == "--csv") csv = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    const Medium in{eps1, 1.0, chi1}, out{1.0, 1.0, 0.0}; const cplx k2 = om;
    std::ofstream f; if (!csv.empty()) { f.open(csv, std::ios::app); f.seekp(0, std::ios::end); if (f.tellp() == 0) f << "n,N,unknowns,omega,eps1_re,eps1_im,chi_re,chi_im,Qplus,Qminus,CD,its_plus,its_minus,MB_H,t_build_s,t_solve_s\n"; }
    std::printf("%4s %7s %8s %11s %11s %11s %9s %7s %9s %8s\n", "n", "N", "Unbek.", "Q_+", "Q_-", "CD/(pi a^2)", "It. +/-", "MB H", "Aufbau s", "Loesen s");
    for (int n : ns) {
        TriangleMesh m = make_icosphere(n);
        auto t0 = std::chrono::steady_clock::now();
        HMatrixParams p; p.eps = heps;
        KernelEntries Ep(m, in.k(om, +1)), Em(m, in.k(om, -1)), Eo(m, k2);
        KernelHMatrix Hp(Ep, p), Hm(Em, p), Ho(Eo, p);
        CauchyOperator Cp(m, Hp), Cm(m, Hm), Co(m, Ho); ChiralCauchyOperator E1(Cp, Cm);
        TransmissionOperator T(m, E1, Co, in, out);
        double tb = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        LinOp A = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.apply(x, y); };
        LinOp M = [&](const std::vector<cplx>& x, std::vector<cplx>& y) { T.precondition(x, y); };
        t0 = std::chrono::steady_clock::now();
        real Q[2]; int its[2]; const Vec3 d(0, 0, 1);
        for (int q = 0; q < 2; ++q) {
            int s = q == 0 ? +1 : -1; CVec3 pol = circular_polarization(d, s);
            auto b = project_plane_wave(m, k2, 1.0, d, pol);
            std::vector<cplx> h; GmresResult r = gmres(A, b, h, &M, tol, 200, 3000); its[q] = r.iterations;
            for (std::size_t i = 0; i < h.size(); ++i) h[i] -= b[i];
            Q[q] = extinction_cross_section(m, h, k2, 1.0, d, pol) / pi;
        }
        double ts = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        double MB = (Hp.stats().bytes() + Hm.stats().bytes() + Ho.stats().bytes()) / 1048576.0;
        std::printf("%4d %7zu %8zu %11.6f %11.6f %11.6f %4d/%-4d %7.0f %9.1f %8.1f\n", n, m.size(), 8 * m.size(), Q[0], Q[1], Q[0] - Q[1], its[0], its[1], MB, tb, ts);
        std::fflush(stdout);
        if (f) { f << n << ',' << m.size() << ',' << 8 * m.size() << ',' << om << ',' << eps1.real() << ',' << eps1.imag() << ',' << chi1.real() << ',' << chi1.imag() << ','
                   << Q[0] << ',' << Q[1] << ',' << Q[0] - Q[1] << ',' << its[0] << ',' << its[1] << ',' << MB << ',' << tb << ',' << ts << '\n'; f.flush(); }
    }
    return 0;
}
