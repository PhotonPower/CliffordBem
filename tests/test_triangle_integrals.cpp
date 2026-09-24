#include "cbem/kernel/triangle_integrals.hpp"
#include "cbem/geometry/quadrature.hpp"
#include "check.hpp"
using namespace cbem;
int main() {
    std::array<Vec3, 3> p = {Vec3(0.1, 0.0, 0.2), Vec3(1.0, 0.2, 0.1), Vec3(0.3, 0.9, 0.3)};
    Vec3 cr = cross(p[1] - p[0], p[2] - p[0]); real A = 0.5 * norm(cr); Vec3 n = cr / (2 * A);
    QuadRule fine = QuadRule::subdivided(60);
    for (Vec3 x : {Vec3(0.5, 0.4, 1.0), Vec3(1.5, -0.5, -0.3), Vec3(0.4, 0.3, 0.6)}) {
        Vec3 Ig; real Ii; triangle_integrals(x, p, n, Ig, Ii);
        Vec3 g{0, 0, 0}; real iv = 0;
        for (std::size_t q = 0; q < fine.w.size(); ++q) {
            Vec3 y = p[0] * fine.bary[q][0] + p[1] * fine.bary[q][1] + p[2] * fine.bary[q][2];
            Vec3 d = x - y; real r = norm(d); g += d * (fine.w[q] * A / (r * r * r)); iv += fine.w[q] * A / r;
        }
        CHECK(norm(Ig - g) / norm(g) < 1e-10, "I_grad: %.2e", norm(Ig - g) / norm(g));
        CHECK(std::abs(Ii - iv) / iv < 1e-10, "I_inv: %.2e", std::abs(Ii - iv) / iv);
    }
    REPORT();
}
