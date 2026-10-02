// Prototyp: Singularitaetssubtraktion zweiter Ordnung im Innenintegral gekruemmter Elemente.
// Innenintegral fuer einen Punkt x ausserhalb von Element j (3 Gewichte lambda_b x 7 Komponenten):
//   I_b(x) = int lambda_b(u) Phi_k(x - X(u)) N(u) du,  N = X_u x X_v = J n.
// Stufe 1 (v0.49): Abzug von lambda_b(u) Phi_sing(z_a) N*, z_a = x - X_aff(u), analytisch (triangle_integrals_linear).
// Stufe 2 (hier): zusaetzlich lambda_b(u*) [K0(z_a) (N_u du + N_v dv) - DK0(z_a)[q(d)] N*], q = D du^2 + E du dv + F dv^2;
//   W_k = int d_k K0(z_a) du aus triangle_integrals_linear, V = int DK0(z_a)[q] du hier noch per feiner Quadratur.
// Messung: Fehler des Innenintegrals gegen eine Referenz und Zahl der Korrekturpunkte je Kriterium.
// Ergebnis und Deutung: docs/results_curved.md (v0.52). Bauen (aus dem Projektverzeichnis, nach dem Bau von build/):
//   g++ -std=c++17 -O3 -fopenmp -Iinclude prototype/curved/second_order_subtraction.cpp build/libcbem.a -o second_order
//   ./second_order [n]   (Ikosaederkugel mit 20 n^2 Elementen, Voreinstellung 4)
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

#include "cbem/geometry/quadratic_mesh.hpp"
#include "cbem/geometry/quadrature.hpp"
#include "cbem/kernel/dirac_kernel.hpp"
#include "cbem/kernel/triangle_integrals.hpp"

using namespace cbem;
using C7 = std::array<cplx, 7>;
using Acc = std::array<C7, 3>;
using Bary = std::array<real, 3>;
using Tri = std::array<Bary, 3>;

static C7 comps7(const Vec3& z, const Vec3& n, cplx s, cplx v) {
    return C7{v * dot(z, n), s * n.x, s * n.y, s * n.z, v * (z.x * n.y - z.y * n.x), v * (z.x * n.z - z.z * n.x), v * (z.y * n.z - z.z * n.y)};
}

struct Geo {
    Vec3 A, B, C, D, E, F;
    Vec3 X(const Bary& l) const { const real u = l[1], v = l[2]; return A + (B + D * u + E * v) * u + (C + F * v) * v; }
    void frame(const Bary& l, Vec3& Xu, Vec3& Xv) const { const real u = l[1], v = l[2]; Xu = B + D * (2 * u) + E * v; Xv = C + E * u + F * (2 * v); }
    Vec3 N(const Bary& l) const { Vec3 Xu, Xv; frame(l, Xu, Xv); return cross(Xu, Xv); }
};

static Geo make_geo(const QuadraticMesh& m, std::size_t t) {
    const auto V = m.flat.vertices(t); const auto& M = m.mid[t]; Geo g;
    g.A = V[0]; g.B = V[0] * -3.0 - V[1] + M[0] * 4.0; g.C = V[0] * -3.0 - V[2] + M[2] * 4.0;
    g.D = (V[0] + V[1]) * 2.0 - M[0] * 4.0; g.E = (V[0] - M[0] + M[1] - M[2]) * 4.0; g.F = (V[0] + V[2]) * 2.0 - M[2] * 4.0;
    return g;
}

static const QuadRule R7 = QuadRule::dunavant7();
// konische Gauss-Produktregel (Duffy), n x n Punkte, exakt bis Grad 2n - 2
static QuadRule stroud(int n) {
    std::vector<real> x(n), w(n);
    for (int i = 0; i < n; ++i) {
        real t = std::cos(pi * (i + 0.75) / (n + 0.5));
        for (int it = 0; it < 100; ++it) {
            real p0 = 1, p1 = t;
            for (int k = 2; k <= n; ++k) { real p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
            const real dp = n * (t * p1 - p0) / (t * t - 1);
            const real dt = p1 / dp; t -= dt; if (std::abs(dt) < 1e-16) break;
        }
        real p0 = 1, p1 = t;
        for (int k = 2; k <= n; ++k) { real p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
        const real dp = n * (t * p1 - p0) / (t * t - 1);
        x[i] = 0.5 * (1 - t); w[i] = 1.0 / ((1 - t * t) * dp * dp);   // auf [0,1]: Gewicht 2/((1-t^2)P'^2) * 1/2
    }
    QuadRule R;
    for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) {
        const real s1 = x[a], t1 = x[b];
        const real l1 = s1, l2 = t1 * (1 - s1);
        R.bary.push_back({1 - l1 - l2, l1, l2}); R.w.push_back(2 * w[a] * w[b] * (1 - s1));
    }
    return R;
}
static QuadRule RX = QuadRule::dunavant7();
static const Tri kRef = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
static long g_points = 0;

// adaptive Integration ueber das Parametergebiet; Kriterium: Umkreisradius des Bildes < ratio * Abstand zu x
template <class F>
static void adapt(const Geo& g, const Vec3& x, const Tri& tri, real aref, int depth, real ratio, int maxdepth, F&& f, const QuadRule* R = &R7) {
    const Vec3 c0 = g.X(tri[0]), c1 = g.X(tri[1]), c2 = g.X(tri[2]), c = (c0 + c1 + c2) / 3.0;
    const real rho = std::max(norm(c0 - c), std::max(norm(c1 - c), norm(c2 - c)));
    const real d = norm(x - c) - rho;
    if (depth < maxdepth && !(rho < ratio * d)) {
        const real l01 = norm(c1 - c0), l12 = norm(c2 - c1), l20 = norm(c0 - c2);
        const int e = (l01 >= l12 && l01 >= l20) ? 0 : (l12 >= l20 ? 1 : 2);
        const auto& A = tri[e]; const auto& B = tri[(e + 1) % 3]; const auto& Cc = tri[(e + 2) % 3];
        const Bary M = {(A[0] + B[0]) / 2, (A[1] + B[1]) / 2, (A[2] + B[2]) / 2};
        adapt(g, x, {A, M, Cc}, aref / 2, depth + 1, ratio, maxdepth, f, R);
        adapt(g, x, {M, B, Cc}, aref / 2, depth + 1, ratio, maxdepth, f, R);
        return;
    }
    for (std::size_t p = 0; p < R->w.size(); ++p) {
        Bary l; for (int k = 0; k < 3; ++k) l[k] = R->bary[p][0] * tri[0][k] + R->bary[p][1] * tri[1][k] + R->bary[p][2] * tri[2][k];
        f(l, R->w[p] * aref);
    }
}

struct Tangent { Bary lam; Vec3 X0, Xu, Xv, N; real J; Vec3 n; };

static Tangent tangent_at(const QuadraticMesh& m, const Geo& g, std::size_t j, const Vec3& x) {
    const auto v = m.flat.vertices(j); const Vec3& nf = m.flat.normal[j];
    const real A2 = dot(cross(v[1] - v[0], v[2] - v[0]), nf);
    Bary l; for (int k = 0; k < 3; ++k) l[k] = dot(cross(v[(k + 2) % 3] - v[(k + 1) % 3], x - v[(k + 1) % 3]), nf) / A2;
    auto clamp = [](Bary& q) { real s = 0; for (auto& c : q) { c = std::max(c, 0.0); s += c; } for (auto& c : q) c /= s; };
    clamp(l);
    for (int it = 0; it < 2; ++it) {
        Vec3 Xu, Xv; g.frame(l, Xu, Xv); const Vec3 r = x - g.X(l);
        const real a11 = dot(Xu, Xu), a12 = dot(Xu, Xv), a22 = dot(Xv, Xv), b1 = dot(r, Xu), b2 = dot(r, Xv), det = a11 * a22 - a12 * a12;
        const real du = (a22 * b1 - a12 * b2) / det, dv = (a11 * b2 - a12 * b1) / det;
        l = {l[0] - du - dv, l[1] + du, l[2] + dv}; clamp(l);
    }
    Tangent T; T.lam = l; g.frame(l, T.Xu, T.Xv); T.N = cross(T.Xu, T.Xv); T.J = norm(T.N); T.n = T.N / T.J;
    T.X0 = g.X(l) - T.Xu * l[1] - T.Xv * l[2];
    return T;
}

// voller Integrand an l mit Gewicht w (Referenzmass du dv)
static void add_full(Acc& out, const Geo& g, const Vec3& x, cplx k, const Bary& l, real w) {
    const Vec3 Xp = g.X(l), N = g.N(l), z = x - Xp; const KernelValue kv = dirac_kernel_fast(z, k);
    const C7 c = comps7(z, N, kv.s, kv.vcoef);
    for (int b = 0; b < 3; ++b) for (int q = 0; q < 7; ++q) out[b][q] += w * l[b] * c[q];
}

struct Sub2 { Vec3 Nu, Nv; Vec3 Wu, Wv, V; };   // Ableitungen von N am Fusspunkt, analytische Vektoren

// Innenintegral: order 0 = ohne Subtraktion (nur fuer Referenz), 1 = wie v0.49, 2 = mit Kruemmungsterm
static Acc inner(const QuadraticMesh& m, const Geo& g, std::size_t j, const Vec3& x, cplx k, int order, real ratio, real vratio = 0.02) {
    Acc out; for (auto& a : out) a.fill(cplx(0));
    if (order == 0) { adapt(g, x, kRef, 0.5, 0, ratio, 30, [&](const Bary& l, real w) { ++g_points; add_full(out, g, x, k, l, w); }, ratio < 0.05 ? &R7 : &RX); return out; }
    const Tangent T = tangent_at(m, g, j, x);
    const std::array<Vec3, 3> tri = {T.X0, T.X0 + T.Xu, T.X0 + T.Xv};
    std::array<Vec3, 3> Ig; std::array<real, 3> Ii;
    triangle_integrals_linear(x, tri, T.n, Ig, Ii);
    const cplx ik = cplx(0, 1) * k;
    for (int b = 0; b < 3; ++b) { const C7 c = comps7(Ig[b] / (4 * pi), T.n, -ik * Ii[b] / (4 * pi), 1.0); for (int q = 0; q < 7; ++q) out[b][q] += c[q]; }
    // Ig, Ii sind Flaechenintegrale (dS = J* du): int lambda_b Phi_sing(z_a) n* dS = int lambda_b Phi_sing(z_a) N* du
    Sub2 S{};
    const real u0 = T.lam[1], v0 = T.lam[2];
    auto q_of = [&](real du, real dv) { return g.D * (du * du) + g.E * (du * dv) + g.F * (dv * dv); };
    auto DK0 = [](const Vec3& z, const Vec3& q) { const real r = norm(z); const Vec3 zh = z / r; return (q - zh * (3 * dot(zh, q))) / (4 * pi * r * r * r); };
    if (order == 2) {
        S.Nu = cross(g.D * 2.0, T.Xv) + cross(T.Xu, g.E); S.Nv = cross(g.E, T.Xv) + cross(T.Xu, g.F * 2.0);
        Vec3 Isum = Ig[0] + Ig[1] + Ig[2];
        // int d_k K0(z_a) du = (1/J*) int (lambda_k - lambda_k*) (x - y)/(4 pi r^3) dS
        S.Wu = (Ig[1] - Isum * u0) / (4 * pi * T.J); S.Wv = (Ig[2] - Isum * v0) / (4 * pi * T.J);
        // V = int DK0(z_a)[q(d)] du, hier per feiner Quadratur (Referenz fuer den Prototyp)
        Vec3 V(0, 0, 0);
        const long keep = g_points;
        Geo ga{T.X0, T.Xu, T.Xv, Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0)};
        adapt(ga, x, kRef, 0.5, 0, vratio, 40, [&](const Bary& l, real w) {
            const Vec3 za = x - ga.X(l); V += DK0(za, q_of(l[1] - u0, l[2] - v0)) * w; });
        g_points = keep;
        S.V = V;
        const real lb[3] = {T.lam[0], T.lam[1], T.lam[2]};
        const C7 a = comps7(S.Wu, S.Nu, 0.0, 1.0), b2 = comps7(S.Wv, S.Nv, 0.0, 1.0), c = comps7(S.V, T.N, 0.0, 1.0);
        for (int b = 0; b < 3; ++b) for (int q = 0; q < 7; ++q) out[b][q] += lb[b] * (a[q] + b2[q] - c[q]);
    }
    adapt(g, x, kRef, 0.5, 0, ratio, 16, [&](const Bary& l, real w) {
        ++g_points;
        const Vec3 Xp = g.X(l), N = g.N(l), z = x - Xp; const KernelValue kv = dirac_kernel_fast(z, k);
        const C7 full = comps7(z, N, kv.s, kv.vcoef);
        const Vec3 za = x - (T.X0 + T.Xu * l[1] + T.Xv * l[2]); const real ra = norm(za);
        const C7 sing = comps7(za, T.N, -ik / (4 * pi * ra), 1.0 / (4 * pi * ra * ra * ra));
        C7 t2; t2.fill(cplx(0));
        if (order == 2) {
            const real du = l[1] - u0, dv = l[2] - v0;
            const Vec3 K0 = za / (4 * pi * ra * ra * ra);
            const C7 a = comps7(K0, S.Nu * du + S.Nv * dv, 0.0, 1.0), c = comps7(DK0(za, q_of(du, dv)), T.N, 0.0, 1.0);
            for (int q = 0; q < 7; ++q) t2[q] = a[q] - c[q];
        }
        for (int b = 0; b < 3; ++b) {
            const real lb2 = order == 2 ? T.lam[b] : 0.0;
            for (int q = 0; q < 7; ++q) out[b][q] += w * (l[b] * (full[q] - sing[q]) - lb2 * t2[q]);
        }
    }, ratio < 0.05 ? &R7 : &RX);
    return out;
}

static double relerr(const Acc& a, const Acc& r) {
    double n = 0, d = 0;
    for (int b = 0; b < 3; ++b) for (int q = 0; q < 7; ++q) { n += std::norm(r[b][q]); d += std::norm(a[b][q] - r[b][q]); }
    return std::sqrt(d / n);
}

int main(int argc, char** argv) {
    const QuadraticMesh m = quadratic_icosphere(argc > 1 ? std::atoi(argv[1]) : 4);
    const cplx k(1.3, 0.05);
    const real near_factor = 2.5;   // wird unten aus der Paarliste bestimmt (wie is_near: Schwerpunktabstand)
    // Paare: nahe, nicht benachbart (wie test_curved, i += 7)
    std::vector<std::pair<std::size_t, std::size_t>> P;
    for (std::size_t i = 0; i < m.size(); i += 7)
        for (std::size_t j = 0; j < m.size(); ++j) {
            int c = 0; for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) c += (m.flat.T[i][a] == m.flat.T[j][b]);
            const real d = norm(m.flat.centroid[i] - m.flat.centroid[j]);
            if (c == 0 && d < near_factor * std::max(m.flat.hmax[i], m.flat.hmax[j])) P.push_back({i, j});
        }
    // aeussere Punkte: 4 x 7 Punkte je Element i (unterteilte Regel), dazu Punkte nahe der Kante zu j
    const QuadRule Rs = QuadRule::subdivided(2);
    struct Case { std::size_t j; Vec3 x; Acc ref; Geo g; };
    std::vector<Case> cases;
    for (auto [i, j] : P) {
        const Geo gi = make_geo(m, i);
        for (std::size_t p = 0; p < Rs.w.size(); p += 3) cases.push_back({j, gi.X({Rs.bary[p][0], Rs.bary[p][1], Rs.bary[p][2]}), {}, make_geo(m, j)});
    }
    std::printf("%zu Paare, %zu Punkte x\n", P.size(), cases.size());
    auto t0 = std::chrono::steady_clock::now();
    double dref = 0;
    for (std::size_t c = 0; c < cases.size(); ++c) {
        cases[c].ref = inner(m, cases[c].g, cases[c].j, cases[c].x, k, 1, 0.03);
        if (c % 25 == 0) dref = std::max(dref, relerr(inner(m, cases[c].g, cases[c].j, cases[c].x, k, 0, 0.04), cases[c].ref));
    }
    std::printf("Referenz (Stufe 1, 0,03) gegen ungeteilte adaptive Quadratur (0,04): %.1e  [%.1f s]\n", dref,
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    struct Rl { const char* name; QuadRule r; };
    std::vector<Rl> rules = {{"D7", QuadRule::dunavant7()}, {"S4", stroud(4)}, {"S5", stroud(5)}, {"S6", stroud(6)}};
    { double sw = 0; for (real w : rules[2].r.w) sw += w; std::printf("Gewichtssumme S5 %.15f\n", sw); }
    for (auto& rl : rules) {
        RX = rl.r;
        for (int order : {0, 1, 2})
            for (real ratio : {0.5, 1.0, 2.0, 100.0}) {
                g_points = 0; double w = 0, s2 = 0;
                for (auto& cs : cases) { const double e = relerr(inner(m, cs.g, cs.j, cs.x, k, order, ratio), cs.ref); w = std::max(w, e); s2 += e * e; }
                std::printf("%s Stufe %d, Kriterium %5.1f: max %.1e, rms %.1e, %.1f Punkte je x\n", rl.name, order, ratio, w, std::sqrt(s2 / cases.size()),
                            double(g_points) / cases.size());
            }
    }
}
