// Duennschicht-Naeherung erster Ordnung (eine Flaeche): Flaechenoperatoren, Rueckfuehrung auf T_1, neutrale Schicht
// exakt ohne Wirkung, Schichtwirkung gegen Aden-Kerker (tools/mie_coated.py); Dirac-Form 2. Ordnung (Formoperator,
// Fehler O((d/a)^2)).
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
    const auto r0 = ThinLayerScatteringProblem(s, gold, {}, om, {}, 0.0, hp, EntryParams{}, ThinLayerModel::Jump1).solve_plane_wave(d, px, so);
    // 2. ohne Schicht = T_1; Schicht aus Aussenmedium: exakt keine Wirkung
    {
        const auto b = ScatteringProblem({s}, {gold}, om, {}, hp).solve_plane_wave(d, px, so);
        CHECK(std::abs(r0.sigma_ext - b.sigma_ext) < 1e-9 * b.sigma_ext, "ohne Schicht nicht T_1");
        const auto rn = ThinLayerScatteringProblem(s, gold, {Coating{0.05, Medium{}}}, om, {}, 0.0, hp, EntryParams{}, ThinLayerModel::Jump1).solve_plane_wave(d, px, so);
        CHECK(std::abs(rn.sigma_ext - r0.sigma_ext) < 1e-12 * r0.sigma_ext, "neutrale Schicht wirkt");
    }
    // 3. Glasschale d = 0,01 auf dem groben Netz (d/h = 0,04): Wirkung gegen Aden-Kerker (exakt 0,090008 bzw. 0,0040723)
    {
        const auto r = ThinLayerScatteringProblem(s, gold, {Coating{0.01, Medium{2.25, 1.0, 0.0}}}, om, {}, 0.0, hp, EntryParams{}, ThinLayerModel::Jump1).solve_plane_wave(d, px, so);
        const real ds = r.sigma_ext - r0.sigma_ext, dp = std::arg(r.forward) - std::arg(r0.forward);
        std::printf("  d = 0,01, n = 4: Delta sigma %.6f (exakt 0,090008), Delta arg S %.7f (exakt 0,0040723), %d It.\n", ds, dp, r.iterations);
        CHECK(std::abs(ds / 0.090008 - 1) < 0.08, "Delta sigma zu ungenau");
        CHECK(std::abs(dp / 0.0040723 - 1) < 0.08, "Delta arg S zu ungenau");
        // lineare Abhaengigkeit von d (erste Ordnung): halbe Dicke -> halbe Wirkung
        const auto rh = ThinLayerScatteringProblem(s, gold, {Coating{0.005, Medium{2.25, 1.0, 0.0}}}, om, {}, 0.0, hp, EntryParams{}, ThinLayerModel::Jump1).solve_plane_wave(d, px, so);
        const real q = (rh.sigma_ext - r0.sigma_ext) / ds;
        std::printf("  Wirkung(d/2) / Wirkung(d) = %.4f\n", q);
        CHECK(std::abs(q - 0.5) < 0.02, "Wirkung nicht naeherungsweise linear in d");
    }
    // 4. Formoperator (glatte Normalen) auf der Kugel mit Radius 2: S = P / 2, Fehler O(h)
    {
        real e[2]; int k = 0;
        for (int n : {6, 12}) {
            const TriangleMesh m = make_icosphere(n, 2.0); SurfaceFV fv(m); e[k] = 0;
            for (std::size_t t = 0; t < m.size(); ++t) { const Vec3& nn = m.normal[t];
                for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) e[k] = std::max(e[k], std::abs(fv.shape[t][3 * a + b] - ((a == b) - nn[a] * nn[b]) / 2.0)); }
            ++k;
        }
        std::printf("  Formoperator: max. Fehler %.2e (n=6), %.2e (n=12)\n", e[0], e[1]);
        CHECK(e[1] < 0.02 && e[1] < 0.7 * e[0], "Formoperator inkonsistent");
    }
    // 5. Dirac-Form 2. Ordnung: ohne Schicht T_1, neutrale Schicht ohne Wirkung, Schichtwirkung gegen Aden-Kerker
    {
        const auto P = [&](const std::vector<Coating>& cs) { return ThinLayerScatteringProblem(s, gold, cs, om, {}, 0.0, hp, EntryParams{}, ThinLayerModel::Dirac2).solve_plane_wave(d, px, so); };
        const auto z = P({});
        CHECK(std::abs(z.sigma_ext - r0.sigma_ext) < 1e-12 * r0.sigma_ext, "Dirac2 ohne Schicht nicht T_1");
        const auto nn = P({Coating{0.05, Medium{}}});
        CHECK(std::abs(nn.sigma_ext - r0.sigma_ext) < 1e-12 * r0.sigma_ext, "Dirac2: neutrale Schicht wirkt");
        const real exs[2] = {0.483026, 1.053497}, exp_[2] = {0.0208498, 0.0428810}, dd[2] = {0.05, 0.1}, tol[2] = {0.02, 0.05};
        for (int i = 0; i < 2; ++i) {
            const auto r = P({Coating{dd[i], Medium{2.25, 1.0, 0.0}}});
            const real ds = r.sigma_ext - r0.sigma_ext, dp = std::arg(r.forward) - std::arg(r0.forward);
            std::printf("  Dirac2 d = %.2f, n = 4: Delta sigma %+.2f %%, Delta arg S %+.2f %% gegen Aden-Kerker, %d It.\n", dd[i], 100 * (ds / exs[i] - 1), 100 * (dp / exp_[i] - 1), r.iterations);
            CHECK(std::abs(ds / exs[i] - 1) < tol[i] && std::abs(dp / exp_[i] - 1) < tol[i], "Dirac2 zu ungenau bei d = %.2f", dd[i]);
        }
    }
    // 6. Dirac2Fit (quadratische Anpassung): Laplace-Beltrami punktweise konsistent, Schicht nach innen (Extrapolation im Gold)
    {
        real e[2]; int k = 0;
        for (int n : {6, 12}) {
            const TriangleMesh m = make_icosphere(n); SurfaceFV fv(m); std::vector<cplx> phi(m.size()), lap; e[k] = 0;
            for (std::size_t t = 0; t < m.size(); ++t) { const Vec3 c = m.centroid[t] / norm(m.centroid[t]); phi[t] = 3 * c.z * c.z - 1; }
            fv.laplace_fit(phi, lap);
            for (std::size_t t = 0; t < m.size(); ++t) e[k] = std::max(e[k], std::abs(lap[t] + 6.0 * phi[t]));
            ++k;
        }
        std::printf("  Laplace-Beltrami (quadratische Anpassung): max. Fehler %.3f (n=6), %.3f (n=12)\n", e[0], e[1]);
        CHECK(e[1] < 0.06 && e[1] < 0.5 * e[0], "Laplace-Beltrami nicht punktweise konsistent");
        const TriangleMesh s8 = make_icosphere(8);
        const auto b8 = ThinLayerScatteringProblem(s8, gold, {}, om, {}, 0.0, hp).solve_plane_wave(d, px, so);
        const auto r8 = ThinLayerScatteringProblem(s8, gold, {Coating{0.02, Medium{}}}, om, {}, 1.0, hp, EntryParams{}, ThinLayerModel::Dirac2Fit).solve_plane_wave(d, px, so);
        const real ds = r8.sigma_ext - b8.sigma_ext;             // Vakuumschicht nach innen = Kugel mit Radius 0,98 (exakt -0,208992)
        std::printf("  Dirac2Fit, Vakuum nach innen d = 0,02, n = 8: Delta sigma %+.2f %% gegen exakt\n", 100 * (ds / -0.208992 - 1));
        CHECK(std::abs(ds / -0.208992 - 1) < 0.01, "Dirac2Fit: Schicht nach innen zu ungenau");
    }
    // 7. mehrere Koerper: weit getrennt additiv
    {
        const TriangleMesh a = translated(s, Vec3(-30, 0, 0)), b = translated(s, Vec3(30, 0, 0));
        const Medium glass{2.25, 1.0, 0.0};
        const auto two = ThinLayerScatteringProblem({ThinBody{a, gold, {Coating{0.05, glass}}}, ThinBody{b, gold, {Coating{0.05, glass}}}}, om, {}, hp).solve_plane_wave(d, px, so);
        const auto one = ThinLayerScatteringProblem(s, gold, {Coating{0.05, glass}}, om, {}, 0.0, hp).solve_plane_wave(d, px, so);
        std::printf("  zwei Koerper im Abstand 60: sigma %.6f, 2 x einzeln %.6f\n", two.sigma_ext, 2 * one.sigma_ext);
        CHECK(std::abs(two.sigma_ext / (2 * one.sigma_ext) - 1) < 0.02, "mehrere Koerper nicht additiv");
    }
    REPORT();
}
