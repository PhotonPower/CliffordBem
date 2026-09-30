#include "cbem/geometry/mesh.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>

namespace cbem {

void TriangleMesh::compute_geometry() {
    const std::size_t N = T.size();
    normal.resize(N); centroid.resize(N); area.resize(N); hmax.resize(N);
    for (std::size_t t = 0; t < N; ++t) {
        auto v = vertices(t);
        Vec3 cr = cross(v[1] - v[0], v[2] - v[0]);
        real a2 = norm(cr);
        area[t] = 0.5 * a2; normal[t] = cr / a2;
        centroid[t] = (v[0] + v[1] + v[2]) / 3.0;
        hmax[t] = std::max({norm(v[1] - v[0]), norm(v[2] - v[1]), norm(v[0] - v[2])});
    }
}

void TriangleMesh::orient_outward(const Vec3& center) {
    for (auto& tri : T) {
        Vec3 a = P[tri[0]], b = P[tri[1]], c = P[tri[2]];
        Vec3 n = cross(b - a, c - a); Vec3 m = (a + b + c) / 3.0 - center;
        if (dot(n, m) < 0) std::swap(tri[1], tri[2]);
    }
    compute_geometry();
}

namespace {
struct KeyCmp {
    bool operator()(const std::array<long long, 3>& a, const std::array<long long, 3>& b) const { return a < b; }
};
std::array<long long, 3> key_of(const Vec3& p) {
    auto q = [](real v) { return static_cast<long long>(std::llround(v * 1e10)); };
    return {q(p.x), q(p.y), q(p.z)};
}
}  // namespace

TriangleMesh translated(const TriangleMesh& m, const Vec3& shift, real scale) {
    TriangleMesh r = m;
    for (auto& p : r.P) p = p * scale + shift;
    r.compute_geometry();
    return r;
}

MultiBodyMesh make_multibody(const std::vector<TriangleMesh>& parts) {
    MultiBodyMesh mb; mb.parts = parts; mb.body_begin.push_back(0);
    for (const auto& p : parts) {
        const int off = static_cast<int>(mb.all.P.size());
        mb.all.P.insert(mb.all.P.end(), p.P.begin(), p.P.end());
        for (auto t : p.T) mb.all.T.push_back({t[0] + off, t[1] + off, t[2] + off});
        mb.body_begin.push_back(mb.all.T.size());
    }
    mb.all.compute_geometry();
    return mb;
}

TriangleMesh make_icosphere(int n, real radius) {
    const real t = (1.0 + std::sqrt(5.0)) / 2.0;
    std::vector<Vec3> V0 = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t},
                            {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
    for (auto& v : V0) v = v / norm(v);
    const int F0[20][3] = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4},
                           {11, 10, 2}, {10, 7, 6}, {7, 1, 8}, {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8},
                           {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
    TriangleMesh m; std::map<std::array<long long, 3>, int, KeyCmp> ids;
    auto vid = [&](Vec3 p) {
        p = p * (radius / norm(p)); auto k = key_of(p);
        auto it = ids.find(k); if (it != ids.end()) return it->second;
        int id = static_cast<int>(m.P.size()); m.P.push_back(p); ids[k] = id; return id;
    };
    for (auto& f : F0) {
        Vec3 A = V0[f[0]], B = V0[f[1]], C = V0[f[2]];
        std::map<std::pair<int, int>, int> idx;
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n - i; ++j)
                idx[{i, j}] = vid(A + (B - A) * (real(i) / n) + (C - A) * (real(j) / n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n - i; ++j) {
                m.T.push_back({idx[{i, j}], idx[{i + 1, j}], idx[{i, j + 1}]});
                if (j < n - i - 1) m.T.push_back({idx[{i + 1, j}], idx[{i + 1, j + 1}], idx[{i, j + 1}]});
            }
    }
    m.orient_outward();
    return m;
}

namespace {
TriangleMesh cube_from_nodes(const std::vector<real>& s);
}

TriangleMesh make_cube_uniform(int n) {
    std::vector<real> s(n + 1);
    for (int i = 0; i <= n; ++i) s[i] = -1.0 + 2.0 * i / n;
    return cube_from_nodes(s);
}

TriangleMesh make_cube_graded(int L) {
    std::vector<real> half = {0.0, 0.5};
    for (int k = 2; k <= L; ++k) half.push_back(1.0 - std::pow(0.5, k));
    half.push_back(1.0);
    std::sort(half.begin(), half.end()); half.erase(std::unique(half.begin(), half.end()), half.end());
    std::vector<real> s;
    for (auto it = half.rbegin(); it != half.rend(); ++it) s.push_back(-*it);
    for (std::size_t i = 1; i < half.size(); ++i) s.push_back(half[i]);
    return cube_from_nodes(s);
}

namespace {
TriangleMesh cube_from_nodes(const std::vector<real>& s) {
    TriangleMesh m; std::map<std::array<long long, 3>, int, KeyCmp> ids;
    auto vid = [&](const Vec3& p) {
        auto k = key_of(p); auto it = ids.find(k); if (it != ids.end()) return it->second;
        int id = static_cast<int>(m.P.size()); m.P.push_back(p); ids[k] = id; return id;
    };
    for (int ax = 0; ax < 3; ++ax)
        for (int sg : {1, -1}) {
            int ua = (ax + 1) % 3, va = (ax + 2) % 3;
            auto pt = [&](real u, real v) { Vec3 p; p[ax] = sg; p[ua] = u; p[va] = v; return vid(p); };
            for (std::size_t i = 0; i + 1 < s.size(); ++i)
                for (std::size_t j = 0; j + 1 < s.size(); ++j) {
                    real uc = 0.5 * (s[i] + s[i + 1]), vc = 0.5 * (s[j] + s[j + 1]);
                    int A = pt(s[i], s[j]), B = pt(s[i + 1], s[j]), C = pt(s[i + 1], s[j + 1]), D = pt(s[i], s[j + 1]);
                    if (uc * vc > 0) { m.T.push_back({A, B, C}); m.T.push_back({A, C, D}); }
                    else { m.T.push_back({A, B, D}); m.T.push_back({B, C, D}); }
                }
        }
    m.orient_outward();
    return m;
}
}  // namespace

namespace {
// Eigenzerlegung einer symmetrischen 3x3-Matrix (Jacobi); A wird zerstoert, V spaltenweise Eigenvektoren
void sym_eigen3(real A[3][3], real lam[3], real V[3][3]) {
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) V[i][j] = (i == j);
    for (int sweep = 0; sweep < 50; ++sweep) {
        real off = std::abs(A[0][1]) + std::abs(A[0][2]) + std::abs(A[1][2]);
        if (off < 1e-15 * (std::abs(A[0][0]) + std::abs(A[1][1]) + std::abs(A[2][2]) + 1e-300)) break;
        for (int p = 0; p < 2; ++p)
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(A[p][q]) < 1e-300) continue;
                real th = 0.5 * (A[q][q] - A[p][p]) / A[p][q];
                real t = (th >= 0 ? 1.0 : -1.0) / (std::abs(th) + std::sqrt(th * th + 1.0));
                real c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k) { real akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq; }
                for (int k = 0; k < 3; ++k) { real apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk; }
                for (int k = 0; k < 3; ++k) { real vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
            }
    }
    for (int i = 0; i < 3; ++i) lam[i] = A[i][i];
}
}  // namespace

namespace {
// Abstand Punkt - Dreieck (naechster Punkt nach Ericson, Real-Time Collision Detection, Abschn. 5.1.5)
real point_tri_dist(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 ab = b - a, ac = c - a, ap = p - a;
    const real d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) return norm(p - a);
    const Vec3 bp = p - b; const real d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return norm(p - b);
    const real vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) return norm(p - (a + ab * (d1 / (d1 - d3))));
    const Vec3 cp = p - c; const real d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return norm(p - c);
    const real vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) return norm(p - (a + ac * (d2 / (d2 - d6))));
    const real va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) return norm(p - (b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)))));
    const real den = 1.0 / (va + vb + vc);
    return norm(p - (a + ab * (vb * den) + ac * (vc * den)));
}
}  // namespace

std::vector<real> distance_to_surface(const TriangleMesh& m, const std::vector<Vec3>& q, real rmax) {
    // Gitter mit Zellgroesse rmax; jedes Dreieck in alle Zellen seiner Huelle (um rmax erweitert) eintragen
    Vec3 lo = m.P[0], hi = m.P[0];
    for (auto& p : m.P) { lo = Vec3(std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)); hi = Vec3(std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)); }
    real hm = 0; for (auto h : m.hmax) hm += h; hm /= std::max<std::size_t>(1, m.hmax.size());
    const real c = std::max(rmax, hm); const real r = rmax; lo = lo - Vec3(c, c, c);   // Zellgroesse >= Suchradius
    auto cell = [&](real v, real o) { return static_cast<long>(std::floor((v - o) / c)); };
    std::map<std::array<long, 3>, std::vector<std::size_t>> grid;
    for (std::size_t t = 0; t < m.T.size(); ++t) {
        Vec3 a = m.P[m.T[t][0]], b = a;
        for (int k = 1; k < 3; ++k) { const Vec3& p = m.P[m.T[t][k]]; a = Vec3(std::min(a.x, p.x), std::min(a.y, p.y), std::min(a.z, p.z)); b = Vec3(std::max(b.x, p.x), std::max(b.y, p.y), std::max(b.z, p.z)); }
        for (long i = cell(a.x - r, lo.x); i <= cell(b.x + r, lo.x); ++i)
            for (long j = cell(a.y - r, lo.y); j <= cell(b.y + r, lo.y); ++j)
                for (long k = cell(a.z - r, lo.z); k <= cell(b.z + r, lo.z); ++k) grid[{i, j, k}].push_back(t);
    }
    std::vector<real> out(q.size(), rmax);
    for (std::size_t i = 0; i < q.size(); ++i) {
        auto it = grid.find({cell(q[i].x, lo.x), cell(q[i].y, lo.y), cell(q[i].z, lo.z)});
        if (it == grid.end()) continue;
        for (std::size_t t : it->second)
            out[i] = std::min(out[i], point_tri_dist(q[i], m.P[m.T[t][0]], m.P[m.T[t][1]], m.P[m.T[t][2]]));
    }
    return out;
}

TriangleMesh make_icosphere_graded(int n, const Vec3& pole, real lambda, real r) {
    TriangleMesh m = make_icosphere(n, 1.0);
    const Vec3 e = pole / norm(pole);
    for (auto& P : m.P) {
        const Vec3 x = P / norm(P); const real xp = dot(x, e); const Vec3 xt = x - e * xp;
        if (1 + xp < 1e-12) { P = x * r; continue; }                       // Gegenpol bleibt
        const Vec3 w = xt * (lambda / (1 + xp)); const real w2 = dot(w, w);
        P = (e * ((1 - w2) / (1 + w2)) + w * (2 / (1 + w2))) * r;
    }
    m.compute_geometry();
    return m;
}

real winding_number(const TriangleMesh& m, const Vec3& x) {
    real omega = 0;
    for (const auto& tr : m.T) {                                        // Raumwinkel je Dreieck (Van Oosterom, Strackee 1983)
        const Vec3 a = m.P[tr[0]] - x, b = m.P[tr[1]] - x, c = m.P[tr[2]] - x;
        const real la = norm(a), lb = norm(b), lc = norm(c);
        const real num = dot(a, cross(b, c));
        const real den = la * lb * lc + dot(a, b) * lc + dot(a, c) * lb + dot(b, c) * la;
        omega += 2 * std::atan2(num, den);
    }
    return omega / (4 * pi);
}

void require_separated(const TriangleMesh& a0, const TriangleMesh& b0, real extra, const std::string& what) {
    TriangleMesh ca, cb;                                                // Geometrie bei Bedarf berechnen
    const TriangleMesh& a = a0.centroid.size() == a0.size() ? a0 : (ca = a0, ca.compute_geometry(), ca);
    const TriangleMesh& b = b0.centroid.size() == b0.size() ? b0 : (cb = b0, cb.compute_geometry(), cb);
    real ha = 0, hb = 0; for (real h : a.hmax) ha += h; for (real h : b.hmax) hb += h;
    const real h = 0.5 * (ha / std::max<std::size_t>(1, a.hmax.size()) + hb / std::max<std::size_t>(1, b.hmax.size()));
    auto box = [](const TriangleMesh& m, Vec3& lo, Vec3& hi) {
        lo = hi = m.P[0];
        for (auto& p : m.P) { lo = Vec3(std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)); hi = Vec3(std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)); }
    };
    Vec3 alo, ahi, blo, bhi; box(a, alo, ahi); box(b, blo, bhi);
    const real tol = extra + 0.05 * h;
    if (alo.x > bhi.x + tol || blo.x > ahi.x + tol || alo.y > bhi.y + tol || blo.y > ahi.y + tol || alo.z > bhi.z + tol || blo.z > ahi.z + tol) return;
    for (int side = 0; side < 2; ++side) {                              // Durchdringung: Knoten der einen Flaeche in der anderen
        const TriangleMesh& p = side ? b : a; const TriangleMesh& q = side ? a : b;
        Vec3 lo, hi; box(q, lo, hi);
        for (const Vec3& v : p.P)
            if (v.x >= lo.x && v.x <= hi.x && v.y >= lo.y && v.y <= hi.y && v.z >= lo.z && v.z <= hi.z && winding_number(q, v) > 0.5)
                throw std::runtime_error(what + ": Flaechen durchdringen sich");
    }
    const real cut = extra + 2 * h;                                     // Beruehrung: kleinster Abstand (Knoten und Schwerpunkte)
    real dmin = cut;
    for (int side = 0; side < 2; ++side) {
        const TriangleMesh& p = side ? b : a; const TriangleMesh& q = side ? a : b;
        std::vector<Vec3> pts = p.P; for (auto& c : p.centroid) pts.push_back(c);
        for (real dd : distance_to_surface(q, pts, cut)) dmin = std::min(dmin, dd);
    }
    if (dmin - extra < 0.05 * h)
        throw std::runtime_error(what + ": Flaechen beruehren sich (Abstand " + std::to_string(dmin - extra) + " < 5 % der Elementgroesse "
                                 + std::to_string(h) + ")");
}

void require_separated_all(const std::vector<const TriangleMesh*>& s, const std::vector<real>& extra, const std::string& what) {
    for (std::size_t i = 0; i < s.size(); ++i)
        for (std::size_t j = i + 1; j < s.size(); ++j)
            require_separated(*s[i], *s[j], (extra.empty() ? 0.0 : extra[i] + extra[j]), what + " (Koerper " + std::to_string(i) + " und " + std::to_string(j) + ")");
}

TriangleMesh offset_surface(const TriangleMesh& m0, real d) {
    TriangleMesh m = m0;
    if (m.normal.size() != m.T.size()) m.compute_geometry();
    const std::size_t nv = m.P.size();
    std::vector<std::array<real, 9>> M(nv, std::array<real, 9>{});
    std::vector<Vec3> b(nv, Vec3{});
    for (std::size_t t = 0; t < m.T.size(); ++t) {
        const Vec3& n = m.normal[t];
        for (int a = 0; a < 3; ++a) {
            const int v = m.T[t][a];
            Vec3 e1 = m.P[m.T[t][(a + 1) % 3]] - m.P[v], e2 = m.P[m.T[t][(a + 2) % 3]] - m.P[v];
            real c = dot(e1, e2) / (norm(e1) * norm(e2));
            real w = std::acos(std::max(-1.0, std::min(1.0, c)));          // Innenwinkel als Gewicht
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) M[v][3 * i + j] += w * n[i] * n[j];
            b[v] += n * w;
        }
    }
    for (std::size_t v = 0; v < nv; ++v) {
        real A[3][3], lam[3], V[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) A[i][j] = M[v][3 * i + j];
        sym_eigen3(A, lam, V);
        const real lmax = std::max({lam[0], lam[1], lam[2]});
        if (!(lmax > 0)) continue;                                       // isolierter Knoten
        Vec3 delta{};
        for (int k = 0; k < 3; ++k) {
            if (lam[k] < 0.02 * lmax) continue;                          // Pseudoinverse: kleine Eigenwerte (glatt) weglassen
            Vec3 e(V[0][k], V[1][k], V[2][k]);
            delta += e * (d * dot(e, b[v]) / lam[k]);
        }
        m.P[v] += delta;
    }
    m.compute_geometry();
    for (std::size_t t = 0; t < m.T.size(); ++t)
        if (!(dot(m.normal[t], m0.normal[t]) > 0.5) || !(m.area[t] > 1e-3 * m0.area[t]))
            throw std::runtime_error("offset_surface: Dreieck " + std::to_string(t) + " klappt um oder entartet (|d| zu gross)");
    // Umstuelpen (Versatz durch das Innere hindurch) aendert das Vorzeichen des eingeschlossenen Volumens
    auto vol = [](const TriangleMesh& q) { real v = 0; for (auto& tr : q.T) v += dot(q.P[tr[0]], cross(q.P[tr[1]], q.P[tr[2]])); return v / 6.0; };
    if (!(vol(m) * vol(m0) > 0)) throw std::runtime_error("offset_surface: Flaeche stuelpt sich um (|d| zu gross)");
    // Faltung oder Durchdringung: Auf einer gueltigen Parallelflaeche haben Knoten und Schwerpunkte ueberall etwa den
    // Abstand |d| von der Originalflaeche (auf Gehrung genau |d|, gekruemmt |d| cos(Winkel der Nachbarnormalen)). Kommt ein
    // Punkt einem anderen Teil der Originalflaeche naeher als 0,8 |d|, hat sich die Flaeche gefaltet (konkave Stelle mit
    // Kruemmungsradius < |d|) oder durchdringt sich (Spalt enger als 2 |d|).
    if (d != 0.0) {
        const real lim = 0.8 * std::abs(d);
        std::vector<Vec3> q = m.P; for (auto& c : m.centroid) q.push_back(c);
        const std::vector<real> dist = distance_to_surface(m0, q, lim);
        for (std::size_t i = 0; i < q.size(); ++i)
            if (dist[i] < lim * (1 - 1e-9))
                throw std::runtime_error("offset_surface: Parallelflaeche faltet oder durchdringt sich (Punkt " + std::to_string(i) + " im Abstand "
                                         + std::to_string(dist[i]) + " < 0,8 |d| von der Originalflaeche)");
    }
    return m;
}

}  // namespace cbem
