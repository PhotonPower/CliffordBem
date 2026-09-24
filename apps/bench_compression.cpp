// Asymptotiktest der H-Matrix-Kompression des Cauchy-Operators (AP 3, Hypothese H3).
// Beispiel:  bench_compression --geometry sphere --n 8,12,16,24 --k 1.5 --eps 1e-4 --mode joint --csv out.csv
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "cbem/operators/cauchy_operator.hpp"

using namespace cbem;

static std::vector<int> parse_list(const std::string& s) {
    std::vector<int> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, ',')) v.push_back(std::stoi(t)); return v;
}

int main(int argc, char** argv) {
    std::string geometry = "sphere", mode = "joint", csv = "";
    std::vector<int> ns = {6, 8, 10};
    double kr = 1.5, ki = 0.0, eps = 1e-4, eta = 1.0, sep = 3.0; int leaf = 32, check_rows = 40;
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--geometry") geometry = nxt(); else if (o == "--n") ns = parse_list(nxt());
        else if (o == "--k") kr = std::stod(nxt()); else if (o == "--ki") ki = std::stod(nxt());
        else if (o == "--eps") eps = std::stod(nxt()); else if (o == "--eta") eta = std::stod(nxt());
        else if (o == "--sep") sep = std::stod(nxt()); else if (o == "--leaf") leaf = std::stoi(nxt());
        else if (o == "--mode") mode = nxt(); else if (o == "--check") check_rows = std::stoi(nxt());
        else if (o == "--csv") csv = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    cplx k(kr, ki);
    std::ofstream out;
    if (!csv.empty()) {
        out.open(csv, std::ios::app);
        out.seekp(0, std::ios::end);
        if (out.tellp() == 0) out << "geometry,n,N,unknowns,mode,eps,eta,leaf,sep,k_re,k_im,n_dense,n_lowrank,MB_dense,MB_lowrank,MB_total,MB_dense_kernel_ref,mean_rank,max_rank,t_build_s,t_matvec_s,rel_err\n";
    }
    std::printf("%-7s %8s %9s %10s %10s %10s %7s %9s %9s %9s\n", "Geom.", "N", "Unbek.", "MB dicht", "MB NR", "MB ges.", "Rang", "Aufbau s", "MatVec s", "Fehler");
    for (int n : ns) {
        TriangleMesh m = geometry == "cube" ? make_cube_graded(n) : make_icosphere(n);
        const std::size_t N = m.size();
        KernelEntries E(m, k);
        HMatrixParams p; p.eps = eps; p.eta = eta; p.leaf = leaf; p.sep_factor = sep;
        p.mode = mode == "comp" ? AcaMode::Componentwise : AcaMode::Joint;
        KernelHMatrix H(E, p); CauchyOperator op(m, H);
        std::mt19937 g(7); std::normal_distribution<real> nd;
        std::vector<cplx> x(8 * N); for (auto& v : x) v = cplx(nd(g), nd(g));
        auto t0 = std::chrono::steady_clock::now(); std::vector<cplx> y; op.apply(x, y);
        double tmv = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        // exakte Referenz fuer zufaellige Zeilen
        std::vector<cplx> Z; op.make_Z(x, Z);
        std::uniform_int_distribution<std::size_t> ud(0, N - 1);
        real num = 0, den = 0;
        for (int r = 0; r < check_rows; ++r) {
            std::size_t i = ud(g); cplx yi[8] = {};
            for (std::size_t j = 0; j < N; ++j) { KernelComp K = E.exact(i, j); for (int c = 0; c < 4; ++c) for (int q = 0; q < 8; ++q) yi[q] += K[c] * Z[32 * j + 8 * c + q]; }
            for (int q = 0; q < 8; ++q) { yi[q] *= -2.0 / std::sqrt(m.area[i]); num += std::norm(y[8 * i + q] - yi[q]); den += std::norm(yi[q]); }
        }
        const HStats& s = H.stats();
        double MBd = 16.0 * s.entries_dense / 1048576.0, MBl = 16.0 * s.entries_lowrank / 1048576.0;
        double MBref = 16.0 * 4.0 * double(N) * double(N) / 1048576.0;
        double err = std::sqrt(num / den);
        std::printf("%-7s %8zu %9zu %10.1f %10.1f %10.1f %7.1f %9.1f %9.3f %9.1e\n", geometry.c_str(), N, 8 * N, MBd, MBl, MBd + MBl, s.mean_rank, s.seconds, tmv, err);
        std::fflush(stdout);
        if (out) out << geometry << ',' << n << ',' << N << ',' << 8 * N << ',' << mode << ',' << eps << ',' << eta << ',' << leaf << ',' << sep << ','
                     << kr << ',' << ki << ',' << s.n_dense << ',' << s.n_lowrank << ',' << MBd << ',' << MBl << ',' << MBd + MBl << ',' << MBref << ','
                     << s.mean_rank << ',' << s.max_rank << ',' << s.seconds << ',' << tmv << ',' << err << '\n';
    }
    return 0;
}
