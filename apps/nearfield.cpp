// Nahfeldkarten (v0.24): Feldverstaerkung |E|^2/|E_0|^2 und optische Chiralitaet C/C_0 in einer Ebene ausserhalb der Koerper.
// Geometrie wie in spectrum: --sphere n [--sphere-dimer g] oder --mesh datei.msh, Laengen in nm mit --unit (Radius der Kugel).
// Schichten: --coating "d:Material[:chi];..." (nach aussen) mit --twoport (Standard) oder --thin 0; --host-chi fuer ein chirales
// Aussenmedium (dann --pol circ). Punkte innerhalb eines Koerpers bzw. seiner Huelle werden markiert (inside = 1).
// Beispiel (Gold-Dimer mit chiraler Schicht, Schnitt durch den Spalt):
//   nearfield --sphere 8 --sphere-dimer 4 --unit 20 --materials Au --nbg 1.33 --coating "1:2.25,0:0.01" --lambda 580 \
//             --pol circ --plane xz --extent "-50:50:101,-30:30:61" --csv results/nf_dimer.csv
#include <chrono>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include "cbem/core/materials.hpp"
#include "cbem/geometry/gmsh_io.hpp"
#include "cbem/problems/twoport_layer_problem.hpp"
#include "cbem/sources/near_field.hpp"
using namespace cbem;
static std::vector<std::string> split(const std::string& s, char c) { std::vector<std::string> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, c)) v.push_back(t); return v; }
int main(int argc, char** argv) {
    std::string mesh, mats = "Au", coating, pol = "circ", plane = "xz", extent = "-40:40:81,-40:40:81", csv, datadir = "data/materials";
    int sph = 0; double unit = 20, nbg = 1.0, host_chi = 0, lambda = 530, at = 0, gap = -1, thin = -1, heps = 1e-6, tol = 1e-8; Vec3 dir(0, 0, 1);
    for (int a = 1; a < argc; ++a) {
        std::string o = argv[a]; auto nxt = [&]() { return std::string(argv[++a]); };
        if (o == "--sphere") sph = std::stoi(nxt()); else if (o == "--sphere-dimer") gap = std::stod(nxt()); else if (o == "--mesh") mesh = nxt();
        else if (o == "--unit") unit = std::stod(nxt()); else if (o == "--materials") mats = nxt(); else if (o == "--nbg") nbg = std::stod(nxt());
        else if (o == "--host-chi") host_chi = std::stod(nxt()); else if (o == "--coating") coating = nxt(); else if (o == "--thin") thin = std::stod(nxt());
        else if (o == "--twoport") thin = -1; else if (o == "--lambda") lambda = std::stod(nxt()); else if (o == "--pol") pol = nxt();
        else if (o == "--dir") { auto v = split(nxt(), ','); dir = Vec3(std::stod(v[0]), std::stod(v[1]), std::stod(v[2])); dir = dir / norm(dir); }
        else if (o == "--plane") plane = nxt(); else if (o == "--extent") extent = nxt(); else if (o == "--at") at = std::stod(nxt());
        else if (o == "--heps") heps = std::stod(nxt()); else if (o == "--tol") tol = std::stod(nxt()); else if (o == "--csv") csv = nxt();
        else if (o == "--data") datadir = nxt();
        else { std::printf("unbekannte Option %s\n", o.c_str()); return 1; }
    }
    std::vector<TriangleMesh> parts;
    if (sph > 0 && gap < 0) parts.push_back(make_icosphere(sph));
    else if (sph > 0) { const real sx = 1.0 + 0.5 * gap / unit;
        parts.push_back(translated(make_icosphere(sph), Vec3(-sx, 0, 0))); parts.push_back(translated(make_icosphere(sph), Vec3(sx, 0, 0))); }
    else for (auto& b : read_gmsh(mesh)) parts.push_back(b.mesh);
    const real om = 2 * pi * unit / lambda;
    auto mat = make_material(mats, datadir);
    const Medium core{mat->eps(lambda), 1.0, 0.0}, bg{nbg * nbg, 1.0, host_chi};
    std::vector<Coating> cs;
    for (auto& c : split(coating, ';')) {
        if (c.empty()) continue; auto q = split(c, ':');
        cs.push_back(Coating{std::stod(q[0]) / unit, Medium{make_material(q[1], datadir)->eps(lambda), 1.0, q.size() > 2 ? std::stod(q[2]) : 0.0}});
    }
    const CVec3 p = pol == "circ" ? circular_polarization(dir, +1) : pol == "circ-" ? circular_polarization(dir, -1) : CVec3{1.0, 0.0, 0.0};
    HMatrixParams hp; hp.eps = heps; SolveOptions so; so.tol = tol;
    auto t0 = std::chrono::steady_clock::now();
    TriangleMesh outer, hull; std::vector<cplx> h; int its = 0;          // hull: Huelle fuer die Markierung (Duennschicht)
    if (cs.empty()) {
        ScatteringProblem P(parts, std::vector<Medium>(parts.size(), core), om, bg, hp);
        auto r = P.solve_plane_wave(dir, p, so); outer = P.mesh(); h = r.h; its = r.iterations;
    } else if (thin >= 0) {
        std::vector<ThinBody> tb; for (auto& pm : parts) tb.push_back(ThinBody{pm, core, cs, thin});
        ThinLayerScatteringProblem P(tb, om, bg, hp); auto r = P.solve_plane_wave(dir, p, so);
        outer = P.mesh(); h = r.h; its = r.iterations;
        // Aussenraum beginnt an der Huelle: fuer die Markierung innerer Punkte die Huelle als Parallelflaeche
        real T = 0; for (auto& c : cs) T += c.thickness;
        std::vector<TriangleMesh> hs_; for (auto& pm : parts) hs_.push_back(offset_surface(pm, (1 - thin) * T));
        hull = make_multibody(hs_).all;
    } else {
        std::vector<TwoPortBody> tb; for (auto& pm : parts) tb.push_back(TwoPortBody{TwoPortLayerProblem::layer_surfaces(pm, cs), core, cs});
        TwoPortLayerProblem P(tb, om, bg, hp); auto r = P.solve_plane_wave(dir, p, so);
        outer = P.outer_mesh(); h.assign(r.h.begin(), r.h.begin() + 8 * outer.size()); its = r.iterations;
    }
    const double ts = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    // Gitter in der Ebene (nm), dritte Koordinate --at
    auto e = split(extent, ','); auto A = split(e[0], ':'), B = split(e[1], ':');
    const real a0 = std::stod(A[0]), a1 = std::stod(A[1]), b0 = std::stod(B[0]), b1 = std::stod(B[1]);
    const int na = std::stoi(A[2]), nb = std::stoi(B[2]);
    std::vector<Vec3> pts;
    for (int j = 0; j < nb; ++j) for (int i = 0; i < na; ++i) {
        const real u = (a0 + (a1 - a0) * i / std::max(1, na - 1)) / unit, v = (b0 + (b1 - b0) * j / std::max(1, nb - 1)) / unit, w = at / unit;
        pts.push_back(plane == "xy" ? Vec3(u, v, w) : plane == "yz" ? Vec3(w, u, v) : Vec3(u, w, v));
    }
    t0 = std::chrono::steady_clock::now();
    auto f = exterior_near_field(outer, h, bg, om, dir, p, pts);
    if (!cs.empty() && thin >= 0) for (std::size_t i = 0; i < pts.size(); ++i) f[i].inside = winding_number(hull, pts[i]) > 0.5;
    const double tf = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    real emax = 0, cmin = 1e300, cmax = -1e300;
    std::size_t nclose = 0; for (auto& q : f) nclose += (!q.inside && q.too_close);
    for (auto& q : f) if (!q.inside && !q.too_close) { emax = std::max(emax, q.enhancement); cmin = std::min(cmin, q.chirality); cmax = std::max(cmax, q.chirality); }
    std::printf("lambda %.1f nm, %zu Punkte (%zu zu nah am Netz): Loesen %.1f s (%d It.), Nahfeld %.1f s; max |E|^2/|E0|^2 = %.2f, C/C0 in [%.3f, %.3f]\n",
                lambda, pts.size(), nclose, ts, its, tf, emax, cmin, cmax);
    if (!csv.empty()) {
        std::ofstream o(csv); o.precision(8);
        o << "x_nm,y_nm,z_nm,inside,too_close,enhancement,chirality,Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im\n";
        for (std::size_t i = 0; i < pts.size(); ++i)
            o << pts[i].x * unit << ',' << pts[i].y * unit << ',' << pts[i].z * unit << ',' << int(f[i].inside) << ',' << int(f[i].too_close) << ',' << f[i].enhancement << ',' << f[i].chirality
              << ',' << f[i].E[0].real() << ',' << f[i].E[0].imag() << ',' << f[i].E[1].real() << ',' << f[i].E[1].imag() << ',' << f[i].E[2].real() << ',' << f[i].E[2].imag() << '\n';
    }
    return 0;
}
