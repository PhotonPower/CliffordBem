// Materialtabellen (Johnson-Christy) und Lebedev-Richtungsquadraturen.
#include "cbem/core/materials.hpp"
#include "cbem/sources/orientation.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    auto au = TabulatedMaterial::from_yaml(CBEM_DATA_DIR "/Au_Johnson.yml");
    CHECK(std::abs(au->lambda_min() - 187.9) < 1e-9 && au->lambda_max() > 1900, "Tabellenbereich %.1f..%.1f", au->lambda_min(), au->lambda_max());
    cplx e = au->eps(187.9); CHECK(std::abs(e - cplx(1.28, 1.188) * cplx(1.28, 1.188)) < 1e-12, "Tabellenpunkt");
    cplx e2 = au->eps(0.5 * (187.9 + 191.6)); cplx nk(0.5 * (1.28 + 1.32), 0.5 * (1.188 + 1.203));
    CHECK(std::abs(e2 - nk * nk) < 1e-12, "lineare Interpolation in n, k");
    cplx e6 = au->eps(600); std::printf("  Au (Johnson-Christy) bei 600 nm: eps = %.3f%+.3fi\n", e6.real(), e6.imag());
    CHECK(e6.real() < -8 && e6.real() > -11 && e6.imag() > 0.5 && e6.imag() < 2, "Plausibilitaet eps(600 nm)");
    bool thrown = false; try { au->eps(100); } catch (const std::out_of_range&) { thrown = true; } CHECK(thrown, "ausserhalb der Tabelle");
    // Lebedev: Momente (1/4pi) int x^a y^b z^c
    auto mom = [](const std::vector<Direction>& D, int a, int b, int c) { real s = 0; for (auto& d : D) s += d.w * std::pow(d.d.x, a) * std::pow(d.d.y, b) * std::pow(d.d.z, c); return s; };
    for (int np : {6, 14, 26}) {
        auto D = lebedev(np);
        CHECK(std::abs(mom(D, 0, 0, 0) - 1) < 1e-14 && std::abs(mom(D, 2, 0, 0) - 1.0 / 3) < 1e-14 && std::abs(mom(D, 1, 1, 0)) < 1e-14, "Lebedev %d: Grad 2", np);
        if (np >= 14) CHECK(std::abs(mom(D, 4, 0, 0) - 0.2) < 1e-14 && std::abs(mom(D, 2, 2, 0) - 1.0 / 15) < 1e-14, "Lebedev %d: Grad 4", np);
        if (np >= 26) CHECK(std::abs(mom(D, 6, 0, 0) - 1.0 / 7) < 1e-14 && std::abs(mom(D, 2, 2, 2) - 1.0 / 105) < 1e-14, "Lebedev %d: Grad 6", np);
    }
    REPORT();
}
