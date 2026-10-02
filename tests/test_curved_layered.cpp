// Geschichtete Koerper auf gekruemmten Elementen (v0.62): Parallelflaeche quadratischer Netze, Rueckfuehrung auf
// CurvedScatteringProblem ohne Schicht (ein und zwei Koerper), Goldkern mit Glasschale (dick, duenn, Doppelschale) gegen
// Aden-Kerker (tools/mie_coated.py), optisches Theorem, chirale Schale gegen tools/mie_chiral_layered.py und
// Spiegelsymmetrie, Nahquadratur mit Randabstand gegen die Voreinstellung.
#include <cmath>
#include <cstdio>

#include "check.hpp"
#include "cbem/problems/curved_layered_problem.hpp"
#include "cbem/problems/curved_problem.hpp"
#include "cbem/sources/fields.hpp"

using namespace cbem;

int main() {
    const real om = 0.5;
    const Medium gold{cplx(-11, 1.2)}, glass{2.25}, vac{};
    const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0};
    SolveOptions so; so.tol = 1e-10;
    const QuadraticMesh q4 = quadratic_icosphere(4);
    // 1. Parallelflaeche: Kugel -> Radius 1 + d (Knoten), zu grosser Versatz nach innen wirft
    {
        const real off = 0.1;
        const QuadraticMesh o = offset_surface(q4, off);
        real e = 0;
        for (const Vec3& p : o.flat.P) e = std::max(e, std::abs(norm(p) - 1 - off));
        for (const auto& M : o.mid) for (const Vec3& p : M) e = std::max(e, std::abs(norm(p) - 1 - off));
        std::printf("  Parallelflaeche (d = 0,1): groesster Radiusfehler der Knoten %.1e\n", e);
        CHECK(e < 1e-4 * off, "Parallelflaeche der Kugel: Radiusfehler %.1e", e);
        bool thrown = false;
        try { offset_surface(q4, -1.2); } catch (const std::exception&) { thrown = true; }
        CHECK(thrown, "zu grosser Versatz nach innen nicht erkannt");
    }
    // 2. ohne Schicht wie CurvedScatteringProblem (ein und zwei Koerper; gleiche Nahquadratur)
    {
        CurvedLayeredGeometry g; add_body(g, q4, gold);
        CurvedLayeredScatteringProblem L(g, om, curved_hmatrix_params(), EntryParams{}, CurvedNearParams{});
        CurvedScatteringProblem C({q4}, {gold}, om);
        const auto a = L.solve_plane_wave(d, px, so); const auto b = C.solve_plane_wave(d, px, so);
        std::printf("  ein Koerper: sigma %.12f / %.12f\n", a.sigma_ext, b.sigma_ext);
        CHECK(std::abs(a.sigma_ext - b.sigma_ext) < 1e-10 * b.sigma_ext && std::abs(a.forward - b.forward) < 1e-10 * std::abs(b.forward),
              "CurvedScatteringProblem nicht reproduziert");
        const QuadraticMesh qa = translated(quadratic_icosphere(3), Vec3(-1.3, 0, 0)), qb = translated(quadratic_icosphere(3), Vec3(1.3, 0, 0));
        CurvedLayeredGeometry g2; add_body(g2, qa, gold); add_body(g2, qb, glass);
        CurvedLayeredScatteringProblem L2(g2, om, curved_hmatrix_params(), EntryParams{}, CurvedNearParams{});
        CurvedScatteringProblem C2({qa, qb}, {gold, glass}, om);
        const auto c = L2.solve_plane_wave(Vec3(1, 0, 1), CVec3{0.0, 1.0, 0.0}, so);
        const auto e = C2.solve_plane_wave(Vec3(1, 0, 1), CVec3{0.0, 1.0, 0.0}, so);
        std::printf("  zwei Koerper: sigma %.12f / %.12f\n", c.sigma_ext, e.sigma_ext);
        CHECK(std::abs(c.sigma_ext - e.sigma_ext) < 1e-10 * e.sigma_ext, "Mehrkoerper: %.12f statt %.12f", c.sigma_ext, e.sigma_ext);
    }
    // 3. Goldkern mit Glasschale gegen Aden-Kerker: dick (bis 1,2), duenn (bis 1,1), Doppelschale (Glas bis 1,1, eps = 4 bis
    //    1,2); Q_ext bezogen auf pi R^2, S(0). Der Imaginaerteil von S(0) ist ungenauer als Q_ext (ohne die Fehlerausloeschung
    //    des optischen Theorems): homogene Goldkugel 1,1e-3 (320 Elemente), 7e-5 (1280)
    {
        struct Case { const char* name; std::vector<real> r; std::vector<Medium> m; real Q; cplx S; real tolQ; };
        for (const Case& c : {Case{"Schale 0,2", {1.2, 1.0}, {glass, gold}, 0.95996079, cplx(0.08639647, -0.30162585), 5e-5},
                              Case{"Schale 0,1", {1.1, 1.0}, {glass, gold}, 0.76475892, cplx(0.05783489, -0.24524718), 5e-4},
                              Case{"Doppelschale", {1.2, 1.1, 1.0}, {Medium{4.0}, glass, gold}, 1.15794913, cplx(0.10421542, -0.32960075), 1e-3}}) {
            std::vector<QuadraticMesh> s; for (real r : c.r) s.push_back(quadratic_icosphere(4, r));
            CurvedLayeredGeometry g; add_layered_body(g, s, c.m);
            CurvedLayeredScatteringProblem P(g, om);
            const auto r = P.solve_plane_wave(d, px, so);
            const real R = c.r[0], Q = r.sigma_ext / (pi * R * R);
            const real eQ = Q / c.Q - 1, eS = std::abs(r.forward - c.S) / std::abs(c.S);
            std::printf("  %-12s: Q_ext %+.2e, S(0) %.2e gegen Aden-Kerker (%d It.)\n", c.name, eQ, eS, r.iterations);
            CHECK(std::abs(eQ) < c.tolQ && eS < 1.5e-3, "%s: Q_ext %.2e, S(0) %.2e", c.name, eQ, eS);
            CHECK(std::abs(4 * pi * r.forward.real() / (om * om) - r.sigma_ext) < 1e-8 * r.sigma_ext, "optisches Theorem: Re S(0) passt nicht zu sigma_ext");
        }
        // beschichteter Koerper ueber die Parallelflaeche (Schale 0,1)
        CurvedLayeredGeometry g; add_coated_body(g, q4, gold, {Coating{0.1, glass}});
        const auto r = CurvedLayeredScatteringProblem(g, om).solve_plane_wave(d, px, so);
        const real eQ = r.sigma_ext / (pi * 1.21) / 0.76475892 - 1;
        std::printf("  add_coated_body (Schale 0,1): Q_ext %+.2e\n", eQ);
        CHECK(std::abs(eQ) < 5e-4, "beschichteter Koerper: Q_ext %.2e", eQ);
    }
    // 4. chirale Glasschale (chi = 0,1) bis 1,2 gegen tools/mie_chiral_layered.py; sigma_s(chi) = sigma_{-s}(-chi)
    {
        const Medium sp{2.25, 1.0, 0.1}, sm{2.25, 1.0, -0.1};
        const real ref[2] = {4.364642220, 4.295021321};
        CurvedLayeredGeometry g; add_layered_body(g, {quadratic_icosphere(4, 1.2), q4}, {sp, gold});
        CurvedLayeredGeometry gm; add_layered_body(gm, {quadratic_icosphere(4, 1.2), q4}, {sm, gold});
        CurvedLayeredScatteringProblem P(g, om), Pm(gm, om);
        const real a = P.solve_plane_wave(d, circular_polarization(d, +1), so).sigma_ext, b = P.solve_plane_wave(d, circular_polarization(d, -1), so).sigma_ext;
        const real c = Pm.solve_plane_wave(d, circular_polarization(d, -1), so).sigma_ext;
        std::printf("  chirale Schale: sigma_+ %+.2e, sigma_- %+.2e gegen Mie, Spiegelsymmetrie %.1e\n", a / ref[0] - 1, b / ref[1] - 1, std::abs(a - c) / a);
        CHECK(std::abs(a / ref[0] - 1) < 1e-4 && std::abs(b / ref[1] - 1) < 1e-4, "chirale Schale weicht von Mie ab");
        CHECK(std::abs(a - c) < 1e-8 * a, "Spiegelsymmetrie verletzt");
        CHECK(std::abs((a - b) / (ref[0] - ref[1]) - 1) < 0.01, "Zirkulardichroismus der Schale: %.4f statt %.4f", a - b, ref[0] - ref[1]);
    }
    // 5. Nahquadratur ueber eine duenne Schicht (d = 0,05, d/h = 0,16): Randabstand gegen die Voreinstellung, Element 0 der
    //    Schale gegen alle nahen Elemente des Kerns
    {
        const QuadraticMesh m = merge_quadratic({quadratic_icosphere(4, 1.05), q4});
        const std::size_t N = q4.size();
        EntryParams ep; ep.cache_near = false;
        const CurvedKernelEntries E0(m, cplx(0.75), ep), E1(m, cplx(0.75), ep, curved_layered_near_params());
        real e = 0, s = 0; int np = 0;
        for (std::size_t j = N; j < 2 * N; ++j) {
            if (!E0.is_near(0, j)) continue;
            ++np;
            const CurvedBlock A = E0.lambda_near(0, j), B = E1.lambda_near(0, j);
            for (int q = 0; q < 9; ++q) for (int c = 0; c < kCurvedComps; ++c) { e = std::max(e, std::abs(A[q][c] - B[q][c])); s = std::max(s, std::abs(A[q][c])); }
        }
        std::printf("  Nahquadratur mit Randabstand, %d Paare ueber die Schicht: Abweichung %.1e\n", np, e / s);
        CHECK(np > 20 && e < 1e-6 * s, "Randabstand aendert die Eintraege: %.2e", e / s);
    }
    REPORT();
}
