// Nahfeld auf gestreckten Elementen: adaptive Aussenregel gegen feine Referenz; Cache gegen direkte Auswertung.
#include "cbem/assembly/kernel_entries.hpp"
#include <random>
#include "check.hpp"
using namespace cbem;
static real rel(const KernelComp& a, const KernelComp& b) { real n = 0, d = 0; for (int c = 0; c < 4; ++c) { n += std::norm(a[c] - b[c]); d += std::norm(b[c]); } return std::sqrt(n / d); }
int main() {
    TriangleMesh m = make_cube_graded(5);
    const cplx k(1.2, 0.3);
    EntryParams pa; pa.cache_near = false;                               // adaptiv, ohne Cache
    EntryParams pr; pr.cache_near = false; pr.adaptive_outer = false; pr.near_subdivision = 16;   // Referenz
    EntryParams pc;                                                      // Standard: adaptiv + Cache
    KernelEntries Ea(m, k, pa), Er(m, k, pr), Ec(m, k, pc);
    std::mt19937 g(11); std::uniform_int_distribution<std::size_t> ud(0, m.size() - 1);
    real emax = 0, ecache = 0; int tested = 0, adj = 0;
    while (tested < 3000) {
        std::size_t i = ud(g), j = ud(g);
        if (!Ea.is_near(i, j)) continue;
        ++tested;
        KernelComp Kr = Ea.adjacency(i, j) != Adjacency::None ? (++adj, Er.exact(i, j)) : Er.near(i, j);
        emax = std::max(emax, rel(Ea.exact(i, j), Kr));
        ecache = std::max(ecache, rel(Ec.exact(i, j), Ea.exact(i, j)));
    }
    std::printf("  %d Nahpaare (%d benachbart): max. Abw. adaptiv/Referenz %.1e, Cache/direkt %.1e; Nahpaare gesamt %zu\n", tested, adj, emax, ecache, Ec.near_pairs());
    CHECK(emax < 2e-4, "Nahfeldregel zu ungenau (%.1e)", emax);
    // Der Cache spiegelt K_ij nach K_ji; die Abweichung zur direkten Auswertung von K_ji ist die
    // Genauigkeit der Nahfeldregeln selbst (halbanalytisch auf gestreckten Elementen ~1e-4)
    CHECK(ecache < 5e-4, "Cache weicht ab (%.1e)", ecache);
    REPORT();
}
