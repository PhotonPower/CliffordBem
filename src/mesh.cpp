#include "cbem/geometry/mesh.hpp"
#include <algorithm>
#include <cmath>
#include <map>
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

TriangleMesh make_cube_graded(int L) {
    std::vector<real> half = {0.0, 0.5};
    for (int k = 2; k <= L; ++k) half.push_back(1.0 - std::pow(0.5, k));
    half.push_back(1.0);
    std::sort(half.begin(), half.end()); half.erase(std::unique(half.begin(), half.end()), half.end());
    std::vector<real> s;
    for (auto it = half.rbegin(); it != half.rend(); ++it) s.push_back(-*it);
    for (std::size_t i = 1; i < half.size(); ++i) s.push_back(half[i]);
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

}  // namespace cbem
