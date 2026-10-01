// Optische Kraefte ueber den Spannungstensor (v0.33): Strahlungsdruck auf eine Goldkugel in Wasser gegen Mie
// (tools/mie_force.py), Uebereinstimmung der drei Wege (Randspuren, Kugel, Parallelflaeche), verschwindende Querkraefte;
// Dimer: Bindungskraefte entgegengesetzt gleich, Summe = Kraft auf die Kugel um beide, Vorzeichen je Polarisation;
// Dipolnaeherung (Stufe 2): ebene Welle = Mie mit n = 1, kleine Glaskugel vor Gold gegen die volle BEM-Rechnung,
// enantioselektive Kraftdifferenz einer kleinen chiralen Kugel gegen die volle BEM-Rechnung.
#include <cmath>
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/optical_force.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    const Medium gold{cplx(-11, 1.2), 1.0, 0.0}, water{1.7689, 1.0, 0.0};
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10; const Vec3 d(0, 0, 1); const CVec3 px{1.0, 0.0, 0.0}, py{0.0, 1.0, 0.0};
    {
        ScatteringProblem P({make_icosphere(8)}, {gold}, 0.5, water, hp); const auto r = P.solve_plane_wave(d, px, so);
        const real ref = 10.629788036;                                         // Mie: sigma_pr 1/2 eps |E0|^2
        const Vec3 Ft = force_from_traces(P.mesh(), r.h, water, {0, P.mesh().size()})[0];
        const Vec3 Fs = force_on_sphere(P.mesh(), r.h, water, 0.5, d, px, Vec3(0, 0, 0), 1.5);
        const Vec3 Fo = force_on_offset(P.mesh(), r.h, water, 0.5, d, px, P.mesh(), 0.2);
        std::printf("  Goldkugel: F_z Spuren %.4f, Kugel %.4f, Parallelflaeche %.4f (Mie %.4f)\n", Ft.z, Fs.z, Fo.z, ref);
        CHECK(std::abs(Fs.z / ref - 1) < 0.025, "Strahlungsdruck weicht von Mie ab");
        CHECK(std::abs(Ft.z / Fs.z - 1) < 0.005 && std::abs(Fo.z / Fs.z - 1) < 0.005, "Wege zur Kraft stimmen nicht ueberein");
        CHECK(std::abs(Fs.x) + std::abs(Fs.y) < 1e-4 * Fs.z, "Querkraft auf die Kugel");
    }
    {
        const real sx = 1.15; const TriangleMesh A = translated(make_icosphere(6), Vec3(-sx, 0, 0)), B = translated(make_icosphere(6), Vec3(sx, 0, 0));
        ScatteringProblem P({A, B}, {gold, gold}, 0.5, water, hp);
        for (int pol = 0; pol < 2; ++pol) {
            const CVec3 p = pol == 0 ? px : py; const auto r = P.solve_plane_wave(d, p, so);
            const Vec3 F1 = force_on_offset(P.mesh(), r.h, water, 0.5, d, p, A, 0.06), F2 = force_on_offset(P.mesh(), r.h, water, 0.5, d, p, B, 0.06);
            const Vec3 Fs = force_on_sphere(P.mesh(), r.h, water, 0.5, d, p, Vec3(0, 0, 0), 3.8);
            std::printf("  Dimer, Polarisation %s: F1 = (%+.3f, %.3f), F2 = (%+.3f, %.3f), Summe F_z %.4f, Kugel um beide %.4f\n",
                        pol == 0 ? "entlang " : "senkrecht", F1.x, F1.z, F2.x, F2.z, F1.z + F2.z, Fs.z);
            CHECK(std::abs(F1.x + F2.x) < 1e-3 * std::abs(F1.x), "Bindungskraefte nicht entgegengesetzt gleich");
            CHECK(std::abs((F1.z + F2.z) / Fs.z - 1) < 0.005, "Summe der Einzelkraefte weicht ab");
            CHECK(pol == 0 ? F1.x > 0 : F1.x < 0, "Vorzeichen der Bindungskraft falsch");
        }
    }
    // Dipolnaeherung
    {
        const CVec3 pc = circular_polarization(d, +1);
        // ebene Welle ueber einen unsichtbaren Koerper: exakt Mie mit n = 1 (Kugel R = 0,3, eps = 4)
        ScatteringProblem I({translated(make_icosphere(4), Vec3(5, 5, 5))}, {water}, 0.5, water, hp); const auto ri = I.solve_plane_wave(d, px, so);
        const auto gi = fields_with_gradients(I.mesh(), ri.h, water, 0.5, d, px, {Vec3(0.2, -0.1, 0.3)}, 1e-3);
        const Vec3 Fi = dipole_particle_force(gi[0], polarizability_from_mie(cplx(2.461520806985727e-06, -0.0015689215238185255),
                                                                           cplx(7.832778966182338e-11, -8.850298845558161e-06), water.k(0.5)), water, 0.5);
        std::printf("  ebene Welle: Dipol F_z %.10e, Mie (n = 1) 9.2276632922e-05\n", Fi.z);
        CHECK(std::abs(Fi.z / 9.227663292185384e-05 - 1) < 1e-6, "Dipolkraft in der ebenen Welle weicht von Mie ab");
        // chirale Glaskugel (R = 0,08, chi = +-0,2) bei (1,4, 0, 0) vor Gold, zirkular: Kraft und Differenz der Enantiomere
        const DipolePolarizability A{cplx(0.0005210437319997587, 5.606082479217233e-09), cplx(-1.4636535189635014e-05, 1.3738518052308093e-09),
                                     cplx(0.0002963869828156086, 2.341654447388878e-09)};   // tools/mie_polarizability.py
        const real rs = 0.08; const Vec3 c(1.4, 0, 0); Vec3 Fb[2];
        for (int e = 0; e < 2; ++e) {
            const TriangleMesh small = translated(make_icosphere(6, rs), c);
            ScatteringProblem Q({make_icosphere(8), small}, {gold, Medium{2.25, 1.0, e == 0 ? 0.2 : -0.2}}, 0.5, water, hp);
            const auto rq = Q.solve_plane_wave(d, pc, so); Fb[e] = force_on_offset(Q.mesh(), rq.h, water, 0.5, d, pc, small, 0.25 * rs);
        }
        ScatteringProblem G({make_icosphere(8)}, {gold}, 0.5, water, hp); const auto rg = G.solve_plane_wave(d, pc, so);
        const auto g = fields_with_gradients(G.mesh(), rg.h, water, 0.5, d, pc, {c}, 1e-3);
        DipolePolarizability Am = A; Am.Ac = -A.Ac;
        const Vec3 Fp = dipole_particle_force(g[0], A, water, 0.5), dD = Fp - dipole_particle_force(g[0], Am, water, 0.5), dB = Fb[0] - Fb[1];
        std::printf("  chirale Kugel vor Gold: |F_Dipol - F_BEM|/|F| %.1e, Differenz der Enantiomere |dD - dB|/|dB| %.1e\n", norm(Fp - Fb[0]) / norm(Fb[0]), norm(dD - dB) / norm(dB));
        CHECK(norm(Fp - Fb[0]) < 0.03 * norm(Fb[0]), "Dipolnaeherung weicht von der vollen BEM-Rechnung ab");
        CHECK(norm(dD - dB) < 0.05 * norm(dB), "enantioselektive Kraft weicht von der vollen BEM-Rechnung ab");
        // v0.35: die chirale Kraft folgt dem Gradienten der optischen Chiralitaet, F(+) - F(-) = -sqrt(eps mu) Re(A_c) grad Im(E*.H)
        Vec3 gC(0, 0, 0); real gc[3];
        for (int i = 0; i < 3; ++i) { cplx q = 0; for (int j = 0; j < 3; ++j) q += std::conj(g[0].dE[i][j]) * g[0].H[j] + std::conj(g[0].E[j]) * g[0].dH[i][j]; gc[i] = std::imag(q); }
        gC = Vec3(gc[0], gc[1], gc[2]) * (-std::sqrt(std::real(water.eps)) * std::real(A.Ac));
        std::printf("  chirale Kraft gegen -sqrt(eps) Re(A_c) grad Im(E*.H): |dD - Modell|/|dD| %.1e\n", norm(dD - gC) / norm(dD));
        CHECK(norm(dD - gC) < 0.05 * norm(dD), "chirale Kraft folgt nicht dem Gradienten der optischen Chiralitaet");
    }
    // v0.36: chirales Aussenmedium (Minkowski-Tensor mit D = eps E + i chi H, B = mu H - i chi E)
    {
        const Medium host{1.7689, 1.0, 0.05}, hostm{1.7689, 1.0, -0.05};
        auto run = [&](const Medium& core, const Medium& out, int s, Vec3& Fs, Vec3& Fo, Vec3& Ff) {
            ScatteringProblem P({make_icosphere(6)}, {core}, 0.5, out, hp); const CVec3 p = circular_polarization(d, s);
            const auto r = P.solve_plane_wave(d, p, so);
            Fs = force_on_sphere(P.mesh(), r.h, out, 0.5, d, p, Vec3(0, 0, 0), 1.5);
            Fo = force_on_offset(P.mesh(), r.h, out, 0.5, d, p, P.mesh(), 0.2);
            Ff = force_from_far_field(P.mesh(), r.h, out, 0.5, d, p, r.sigma_ext); };
        Vec3 a1, b1, c1, a2, b2, c2, a3, b3, c3;
        run(host, host, +1, a1, b1, c1);                                         // unsichtbar
        std::printf("  chirales Aussenmedium, unsichtbare Kugel: |F| = %.1e\n", norm(a1) + norm(b1) + norm(c1));
        CHECK(norm(a1) + norm(b1) + norm(c1) < 1e-10, "Kraft auf eine Kugel aus dem Aussenmedium");
        run(gold, host, +1, a2, b2, c2); run(gold, host, -1, a3, b3, c3);
        std::printf("  Gold in chiralem Wasser (chi 0,05): s=+1 Kugel %.4f, Parallelflaeche %.4f, Fernfeld %.4f; s=-1 %.4f, %.4f, %.4f\n", a2.z, b2.z, c2.z, a3.z, b3.z, c3.z);
        CHECK(std::abs(a2.z / c2.z - 1) < 0.005 && std::abs(b2.z / c2.z - 1) < 0.005, "Spannungstensor im chiralen Medium weicht von der Impulsbilanz ab (s = +1)");
        CHECK(std::abs(a3.z / c3.z - 1) < 0.005 && std::abs(b3.z / c3.z - 1) < 0.005, "Spannungstensor im chiralen Medium weicht von der Impulsbilanz ab (s = -1)");
        CHECK(a3.z > a2.z, "Helizitaetsabhaengigkeit des Strahlungsdrucks fehlt");
        Vec3 a4, b4, c4; run(gold, hostm, -1, a4, b4, c4);                      // Spiegelsymmetrie F_+(chi) = F_-(-chi)
        std::printf("  Spiegelsymmetrie: F_+(chi) %.8f, F_-(-chi) %.8f\n", a2.z, a4.z);
        CHECK(std::abs(a2.z - a4.z) < 1e-6 * a2.z, "Spiegelsymmetrie verletzt");
    }
    REPORT();
}
