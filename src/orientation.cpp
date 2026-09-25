#include "cbem/sources/orientation.hpp"
#include <cmath>
#include <stdexcept>

namespace cbem {

std::vector<Direction> lebedev(int np) {
    std::vector<Direction> r;
    auto axes = [&](real w) { for (int a = 0; a < 3; ++a) for (int s : {1, -1}) { Vec3 v; v[a] = s; r.push_back({v, w}); } };
    auto diag = [&](real w) { const real c = 1 / std::sqrt(3.0); for (int x : {1, -1}) for (int y : {1, -1}) for (int z : {1, -1}) r.push_back({Vec3(x * c, y * c, z * c), w}); };
    auto edges = [&](real w) { const real c = 1 / std::sqrt(2.0);
        for (int a = 0; a < 3; ++a) for (int s1 : {1, -1}) for (int s2 : {1, -1}) { Vec3 v; v[a] = s1 * c; v[(a + 1) % 3] = s2 * c; r.push_back({v, w}); } };
    if (np == 1) r.push_back({Vec3(0, 0, 1), 1.0});
    else if (np == 6) axes(1.0 / 6);
    else if (np == 14) { axes(1.0 / 15); diag(3.0 / 40); }
    else if (np == 26) { axes(1.0 / 21); edges(4.0 / 105); diag(9.0 / 280); }
    else throw std::invalid_argument("lebedev: 1, 6, 14 oder 26 Punkte");
    return r;
}

}  // namespace cbem
