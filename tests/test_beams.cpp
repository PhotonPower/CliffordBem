// Allgemeine einfallende Felder und Kraefte bei Dipol- und Strahlanregung (v0.37): Strahlen als exakte Maxwell-Loesungen,
// Leistung nach Parseval gegen den Poynting-Fluss, ebene Welle ueber die allgemeine Schnittstelle wie bisher, Impulserhaltung
// bei Dipolanregung (Koerper + Emitter + Abstrahlung = 0), kleines Teilchen im fokussierten Strahl gegen die Dipolnaeherung;
// Strahl im chiralen Medium (v0.38); optische Pinzette gegen die GLMT-Referenz (v0.39).
#include <cmath>
#include "cbem/problems/scattering_problem.hpp"
#include "cbem/sources/dipole.hpp"
#include "cbem/sources/optical_force.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    const Medium water{1.7689, 1.0, 0.0};
    // 1. Strahlen (Laengen in um, lambda 1,064 um)
    {
        const real lam = 1.064, om = 2 * pi / lam;
        const BeamField F = BeamField::focused(water, om, Vec3(0, 0, 0), 1.2, 1.0, CVec3{1.0, 0.0, 0.0});
        const BeamField G = BeamField::gaussian(water, om, Vec3(0, 0, 0), 2 * lam, CVec3{1.0, 0.0, 0.0});
        for (const BeamField* B : {&F, &G}) {
            const real h = 1e-4 * lam; real dmax = 0, cmax = 0;
            for (const Vec3 x : {Vec3(0.1, 0.05, 0.02) * lam, Vec3(-0.3, 0.2, 0.4) * lam}) {
                CVec3 E0, H0, Ep[3], Em[3], Hd; B->eval(x, E0, H0);
                for (int i = 0; i < 3; ++i) { Vec3 e(0, 0, 0); (i == 0 ? e.x : i == 1 ? e.y : e.z) = h; B->eval(x + e, Ep[i], Hd); B->eval(x - e, Em[i], Hd); }
                auto d = [&](int i, int j) { return (Ep[i][j] - Em[i][j]) / (2 * h); };
                const cplx div = d(0, 0) + d(1, 1) + d(2, 2); const CVec3 curl{d(1, 2) - d(2, 1), d(2, 0) - d(0, 2), d(0, 1) - d(1, 0)};
                real e2 = 0, r = 0; for (int a = 0; a < 3; ++a) { e2 += std::norm(E0[a]); r += std::norm(curl[a] - cplx(0, om) * H0[a]); }
                dmax = std::max(dmax, std::abs(div) / (om * std::sqrt(e2))); cmax = std::max(cmax, std::sqrt(r / e2) / om);
            }
            std::printf("  Strahl: |div E| %.1e, |rot E - i omega mu H| %.1e\n", dmax, cmax);
            CHECK(dmax < 1e-6 && cmax < 1e-5, "Strahl erfuellt die Maxwell-Gleichungen nicht");
        }
        const real P = G.poynting_flux(0, 8 * lam, 160);
        std::printf("  Gaussstrahl: Poynting-Fluss %.6f (Parseval 1)\n", P);
        CHECK(std::abs(P - 1) < 1e-4, "Leistung des Gaussstrahls");
    }
    // 2. ebene Welle ueber die allgemeine Schnittstelle: dieselbe Projektion wie project_plane_wave
    {
        const TriangleMesh m = make_icosphere(4); const Vec3 d(0.3, 0.2, 0.93); const CVec3 p = circular_polarization(d / norm(d), +1);
        const auto a = PlaneWaveField(water, 0.5, d, p).project(m, water), b = project_plane_wave(m, water.k(0.5), water.eps, d / norm(d), p);
        real e = 0, s = 0; for (std::size_t i = 0; i < a.size(); ++i) { e += std::norm(a[i] - b[i]); s += std::norm(b[i]); }
        std::printf("  ebene Welle allgemein gegen bisher: %.1e\n", std::sqrt(e / s));
        CHECK(std::sqrt(e / s) < 1e-13, "allgemeine ebene Welle weicht ab");
    }
    HMatrixParams hp; hp.eps = 1e-8; SolveOptions so; so.tol = 1e-10;
    // 3. Impulserhaltung bei Dipolanregung
    {
        const Medium gold{cplx(-11, 1.2), 1.0, 0.0}; const Vec3 r0(1.5, 0, 0); const CVec3 p{1.0, 0.0, 0.0};
        ScatteringProblem P({make_icosphere_graded(6, r0, 0.5)}, {gold}, 0.5, water, hp);
        auto src = std::make_shared<DipoleField>(water, 0.5, r0, p);
        const auto b = project_dipole(P.mesh(), water, 0.5, r0, p); const auto r = P.solve_rhs(b, so);
        const Vec3 Fb = force_on_offset(make_near_field_eval(P.mesh(), r.h, b, water, 0.5, src), water, P.mesh(), 0.2);
        const Vec3 Fe = emitter_force(P.mesh(), r.h, b, water, 0.5, *src, 0.005), Pi = radiated_momentum(P.mesh(), r.h, b, water, 0.5, *src);
        const real rel = norm(Fb + Fe + Pi) / norm(Fb);
        std::printf("  Dipolanregung: Kugel %.5f, Emitter %.5f, Abstrahlung %.5f, Summe relativ %.1e\n", Fb.x, Fe.x, Pi.x, rel);
        CHECK(rel < 2e-3, "Impulserhaltung bei Dipolanregung verletzt");
        CHECK(Fb.x > 0 && Fe.x < 0, "Kugel und Emitter ziehen sich nicht an");
    }
    // 4. kleines Teilchen (Polystyrol R = 0,02 um) im fokussierten Strahl: BEM gegen Dipolnaeherung
    {
        const real R = 0.02, om = 2 * pi * R / 1.064; const Medium ps{2.5281, 1.0, 0.0};
        const DipolePolarizability A{cplx(1.5688851148490328, 0.0005061052452322659), cplx(0.004424785413581378, 4.025714163879892e-09), 0.0};
        ScatteringProblem P({make_icosphere(8)}, {ps}, om, water, hp);
        const Vec3 x(0.15, 0.0, 0.1);
        auto beam = std::make_shared<BeamField>(BeamField::focused(water, om, x * (-1.0 / R), 1.2, 1.0, CVec3{1.0, 0.0, 0.0}));
        const auto b = beam->project(P.mesh(), water); const auto r = P.solve_rhs(b, so);
        const Vec3 Fb = force_on_sphere(make_near_field_eval(P.mesh(), r.h, b, water, om, beam), water, Vec3(0, 0, 0), 1.5);
        NearFieldEval only_beam = [&](const std::vector<Vec3>& pts) { std::vector<NearFieldPoint> f(pts.size()); for (std::size_t i = 0; i < pts.size(); ++i) beam->eval(pts[i], f[i].E, f[i].H); return f; };
        const Vec3 Fd = dipole_particle_force(fields_with_gradients(only_beam, {Vec3(0, 0, 0)}, 1e-3)[0], A, water, om);
        std::printf("  kleines Teilchen im Fokus: BEM (%+.4e, %+.4e), Dipol (%+.4e, %+.4e), rel. %.1e\n", Fb.x, Fb.z, Fd.x, Fd.z, norm(Fb - Fd) / norm(Fd));
        CHECK(norm(Fb - Fd) < 0.015 * norm(Fd), "Kraft im fokussierten Strahl weicht von der Dipolnaeherung ab");
    }
    // 5. v0.38: Strahl im chiralen Medium -- chirale Maxwell-Gleichungen, Leistung, Spiegelsymmetrie der Pinzette
    {
        const real lam = 1.064, om = 2 * pi / lam; const Medium mc{1.7689, 1.0, 0.05};
        const real s2 = 1.0 / std::sqrt(2.0);
        const BeamField B = BeamField::focused(mc, om, Vec3(0, 0, 0), 1.2, 0.5, CVec3{s2, cplx(0, s2), 0.0}, 30, 60);
        const real h = 1e-4 * lam; const Vec3 x = Vec3(-0.3, 0.2, 0.4) * lam;
        CVec3 E0, H0, Ep[3], Em[3], Hp[3], Hm[3]; B.eval(x, E0, H0);
        for (int i = 0; i < 3; ++i) { Vec3 e(0, 0, 0); (i == 0 ? e.x : i == 1 ? e.y : e.z) = h; B.eval(x + e, Ep[i], Hp[i]); B.eval(x - e, Em[i], Hm[i]); }
        auto dE = [&](int i, int j) { return (Ep[i][j] - Em[i][j]) / (2 * h); }; auto dH = [&](int i, int j) { return (Hp[i][j] - Hm[i][j]) / (2 * h); };
        const CVec3 cE{dE(1, 2) - dE(2, 1), dE(2, 0) - dE(0, 2), dE(0, 1) - dE(1, 0)}, cH{dH(1, 2) - dH(2, 1), dH(2, 0) - dH(0, 2), dH(0, 1) - dH(1, 0)};
        const cplx I(0, 1); real e2 = 0, r1 = 0, r2 = 0;
        for (int a = 0; a < 3; ++a) { e2 += std::norm(E0[a]); r1 += std::norm(cE[a] - I * om * (mc.mu * H0[a] - I * mc.chi * E0[a])); r2 += std::norm(cH[a] + I * om * (mc.eps * E0[a] + I * mc.chi * H0[a])); }
        const real Pf = B.poynting_flux(0, 6 * lam, 120);
        std::printf("  chiraler Strahl: rot E %.1e, rot H %.1e, Poynting-Fluss %.5f\n", std::sqrt(r1 / e2) / om, std::sqrt(r2 / e2) / om, Pf);
        CHECK(std::sqrt(r1 / e2) / om < 1e-5 && std::sqrt(r2 / e2) / om < 1e-5, "chiraler Strahl erfuellt die Maxwell-Gleichungen nicht");
        CHECK(std::abs(Pf - 1) < 1e-3, "Leistung des chiralen Strahls");
        // Spiegelsymmetrie: Q_z(s = +1, chi) = Q_z(s = -1, -chi) fuer eine achirale Kugel auf der Achse
        const real R = 0.25, omr = 2 * pi * R / lam; real Q[2];
        for (int k = 0; k < 2; ++k) {
            const Medium host{1.7689, 1.0, k == 0 ? 0.05 : -0.05}, ps{2.5281, 1.0, 0.0};
            ScatteringProblem P({make_icosphere(4)}, {ps}, omr, host, hp);
            auto beam = std::make_shared<BeamField>(BeamField::focused(host, omr, Vec3(0, 0, -0.4), 1.2, 1.0, CVec3{s2, cplx(0, k == 0 ? s2 : -s2), 0.0}, 24, 48));
            const auto b = beam->project(P.mesh(), host); const auto r = P.solve_rhs(b, so);
            Q[k] = force_on_sphere(make_near_field_eval(P.mesh(), r.h, b, host, omr, beam), host, Vec3(0, 0, 0), 1.3).z;
        }
        std::printf("  Spiegelsymmetrie der Pinzette: Q_z(+, chi) %.7f, Q_z(-, -chi) %.7f\n", Q[0], Q[1]);
        CHECK(std::abs(Q[0] - Q[1]) < 1e-3 * std::abs(Q[0]), "Spiegelsymmetrie verletzt");
    }
    // 6. v0.39: optische Pinzette gegen die GLMT-Referenz (tools/glmt.py, derselbe Strahl): Polystyrol R = 0,25 um, NA 1,2,
    //    x-polarisiert, Teilchen 0,5 um hinter dem Fokus; GLMT Q_z = -0,0154336, BEM n = 8 etwa 0,3 % darunter
    {
        const real R = 0.25, lam = 1.064, omr = 2 * pi * R / lam; const Medium ps{2.5281, 1.0, 0.0};
        ScatteringProblem P({make_icosphere(8)}, {ps}, omr, water, hp);
        auto beam = std::make_shared<BeamField>(BeamField::focused(water, omr, Vec3(0, 0, -0.5 / R), 1.2, 1.0, CVec3{1.0, 0.0, 0.0}));
        const auto b = beam->project(P.mesh(), water); const auto r = P.solve_rhs(b, so);
        const real Qz = force_on_sphere(make_near_field_eval(P.mesh(), r.h, b, water, omr, beam), water, Vec3(0, 0, 0), 1.3).z / std::sqrt(1.7689);
        const real glmt = -0.01543359404291259;
        std::printf("  Pinzette gegen GLMT: Q_z BEM %.6f, GLMT %.6f (%.2f %%)\n", Qz, glmt, 100 * (Qz / glmt - 1));
        CHECK(std::abs(Qz / glmt - 1) < 0.01, "Pinzette weicht von der GLMT ab");
    }
    REPORT();
}
