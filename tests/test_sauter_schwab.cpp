// Sauter-Schwab: Gewichtssummen, Exaktheit fuer glatte Integranden, Konvergenz fuer die singulaeren
// Kernintegrale im Vergleich mit der unabhaengigen Methode (analytisches Innenintegral, feine Aussenregel).
#include "cbem/assembly/kernel_entries.hpp"
#include "cbem/geometry/quadrature.hpp"
#include "check.hpp"
using namespace cbem;
static real rel(const KernelComp& a, const KernelComp& b) {
    real n = 0, d = 0; for (int c = 0; c < 4; ++c) { n += std::norm(a[c] - b[c]); d += std::norm(b[c]); } return std::sqrt(n / d);
}
int main() {
    // (a) Gewichtssummen und glatter Integrand
    QuadRule fine = QuadRule::subdivided(8);
    auto f = [](const std::array<real, 2>& x, const std::array<real, 2>& y) { return std::exp(0.3 * x[0] + 0.7 * x[1] - 0.5 * y[0] + 0.2 * y[1]); };
    real ref = 0;
    for (std::size_t a = 0; a < fine.w.size(); ++a) for (std::size_t b = 0; b < fine.w.size(); ++b)
        ref += 0.25 * fine.w[a] * fine.w[b] * f({fine.bary[a][1], fine.bary[a][2]}, {fine.bary[b][1], fine.bary[b][2]});
    for (Adjacency ad : {Adjacency::Vertex, Adjacency::Edge, Adjacency::Coincident}) {
        PairRule R = PairRule::sauter_schwab(ad, 5);
        real sw = 0, sf = 0; for (std::size_t q = 0; q < R.w.size(); ++q) { sw += R.w[q]; sf += R.w[q] * f(R.x[q], R.y[q]); }
        CHECK(std::abs(sw - 0.25) < 1e-13, "Gewichtssumme %.15f", sw);
        CHECK(std::abs(sf - ref) / ref < 1e-9, "glatter Integrand: %.2e", std::abs(sf - ref) / ref);
    }
    // (b) singulaere Paare: SS(Ordnung) gegen analytisches Innenintegral mit feiner Aussenregel
    TriangleMesh m = make_icosphere(3);
    const cplx k(1.2, 0.4);
    std::size_t iv = 0, ie = 0, i0 = 0; bool fv = false, fe = false;
    KernelEntries probe(m, k);
    for (std::size_t j = 1; j < m.size() && !(fv && fe); ++j) {
        Adjacency a = probe.adjacency(i0, j);
        if (a == Adjacency::Vertex && !fv) { iv = j; fv = true; }
        if (a == Adjacency::Edge && !fe) { ie = j; fe = true; }
    }
    struct Case { const char* name; std::size_t j; } cases[3] = {{"gleiches Dreieck", i0}, {"gemeinsame Kante", ie}, {"gemeinsame Ecke", iv}};
    for (auto& cs : cases) {
        EntryParams pr; pr.ss_order = 10; KernelEntries Ref(m, k, pr);
        const Adjacency ad = Ref.adjacency(i0, cs.j);
        KernelComp Kr = Ref.sauter_schwab(i0, cs.j, ad);                 // Referenz: Sauter-Schwab Ordnung 10
        std::printf("  %s (Referenz: Sauter-Schwab Ordnung 10):\n", cs.name);
        real prev = 1e9; bool mono = true;
        for (int sub : {4, 16, 48}) { EntryParams p2; p2.sauter_schwab = false; p2.near_subdivision = sub; KernelEntries E(m, k, p2);
            real e = rel(E.near(i0, cs.j), Kr); mono &= (e < prev); prev = e;
            std::printf("     analytisch innen, Aussenregel %2d^2 x 7 Punkte: Abw. %.2e\n", sub, e); }
        CHECK(mono, "unabhaengige Methode konvergiert nicht gegen Sauter-Schwab");
        real e6 = 0;
        for (int ord : {3, 4, 6, 8}) { EntryParams p3; p3.ss_order = ord; KernelEntries E(m, k, p3);
            real e = rel(E.sauter_schwab(i0, cs.j, ad), Kr); if (ord == 6) e6 = e;
            std::printf("     Sauter-Schwab Ordnung %d: Abw. %.2e\n", ord, e); }
        CHECK(e6 < 1e-6, "Sauter-Schwab konvergiert nicht (Ordnung 6: %.1e)", e6);
    }
    REPORT();
}
