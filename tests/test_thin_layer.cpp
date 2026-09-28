// Duennschicht-Naeherung erster Ordnung (eine Flaeche): Flaechenoperatoren, Rueckfuehrung auf T_1, neutrale Schicht
// exakt ohne Wirkung, Schichtwirkung gegen Aden-Kerker (tools/mie_coated.py).
#include <cmath>
#include "cbem/problems/thin_layer_problem.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    // 1. Finite-Volumen-Operatoren auf der Kugel: grad z = e_z - (e_z.n) n,  div_G (a - (a.n) n) = -2 (a.n) / R
    {
        const TriangleMesh m = make_icosphere(12); SurfaceFV fv(m);
        std::vector<cplx> phi(m.size()); std::vector<std::array<cplx, 3>> v(m.size()), g; std::vector<cplx> dv;
        for (std::size_t t = 0; t < m.size(); ++t) {
            phi[t] = m.centroid[t].z; const Vec3& n = m.normal[t];
            for (int d = 0; d < 3; ++d) v[t][d] = (d == 2 ? 1.0 : 0.0) - n.z * n[d];
        }
        fv.grad(phi, g); fv.div(v, dv);
        real eg = 0, ed = 0;
        for (std::size_t t = 0; t < m.size(); ++t) {
            const Vec3& n = m.normal[t];
            for (int d = 0; d < 3; ++d) eg = std::max(eg, std::abs(g[t][d] - ((d == 2 ? 1.0 : 0.0) - n.z * n[d])));
            ed = std::max(ed, std::abs(dv[t] - (-2.0 * n.z)));
        }
        std::printf("  FV auf der Kugel (n=12): max. Fehler grad %.2e, div %.2e\n", eg, ed);
        CHECK(eg < 0.02 && ed < 0.1, "Flaechenoperatoren inkonsistent");
    }
    const real om = 0.5; const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10;
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0};
    const TriangleMesh s = make_icosphere(4);
    const auto r0 = ThinLayerScatteringProblem(s, gold, {}, om, {}, 0.0, hp).solve_plane_wave(d, px, so);
    // 2. ohne Schicht = T_1; Schicht aus Aussenmedium: exakt keine Wirkung
    {
        const auto b = ScatteringProblem({s}, {gold}, om, {}, hp).solve_plane_wave(d, px, so);
        CHECK(std::abs(r0.sigma_ext - b.sigma_ext) < 1e-9 * b.sigma_ext, "ohne Schicht nicht T_1");
        const auto rn = ThinLayerScatteringProblem(s, gold, {Coating{0.05, Medium{}}}, om, {}, 0.0, hp).solve_plane_wave(d, px, so);
        CHECK(std::abs(rn.sigma_ext - r0.sigma_ext) < 1e-12 * r0.sigma_ext, "neutrale Schicht wirkt");
    }
    // 3. Glasschale d = 0,01 auf dem groben Netz (d/h = 0,04): Wirkung gegen Aden-Kerker (exakt 0,090008 bzw. 0,0040723)
    {
        const auto r = ThinLayerScatteringProblem(s, gold, {Coating{0.01, Medium{2.25, 1.0, 0.0}}}, om, {}, 0.0, hp).solve_plane_wave(d, px, so);
        const real ds = r.sigma_ext - r0.sigma_ext, dp = std::arg(r.forward) - std::arg(r0.forward);
        std::printf("  d = 0,01, n = 4: Delta sigma %.6f (exakt 0,090008), Delta arg S %.7f (exakt 0,0040723), %d It.\n", ds, dp, r.iterations);
        CHECK(std::abs(ds / 0.090008 - 1) < 0.08, "Delta sigma zu ungenau");
        CHECK(std::abs(dp / 0.0040723 - 1) < 0.08, "Delta arg S zu ungenau");
        // lineare Abhaengigkeit von d (erste Ordnung): halbe Dicke -> halbe Wirkung
        const auto rh = ThinLayerScatteringProblem(s, gold, {Coating{0.005, Medium{2.25, 1.0, 0.0}}}, om, {}, 0.0, hp).solve_plane_wave(d, px, so);
        const real q = (rh.sigma_ext - r0.sigma_ext) / ds;
        std::printf("  Wirkung(d/2) / Wirkung(d) = %.4f\n", q);
        CHECK(std::abs(q - 0.5) < 0.02, "Wirkung nicht naeherungsweise linear in d");
    }
    REPORT();
}
