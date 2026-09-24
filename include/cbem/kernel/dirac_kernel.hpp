#pragma once
// Fundamentalloesung des Dirac-Operators (nabla - ik): Phi_k(z) = -(nabla + ik) e^{ik|z|}/(4 pi |z|)
//   = Phi(|z|) [ (1/r - ik) z/r - ik ]   (Vektor- plus Skalaranteil).
// Zerlegung: Phi_k = Phi_0 + (-ik/(4 pi r)) + Rem_k mit Phi_0 = z/(4 pi r^3) und beschraenktem Rem_k.
#include "cbem/core/types.hpp"

namespace cbem {

struct KernelValue { cplx s; cplx vcoef; };   // Phi = vcoef * z + s  (vcoef multipliziert den Vektor z)

inline KernelValue dirac_kernel_full(const Vec3& z, cplx k) {
    real r = norm(z);
    cplx Phi = std::exp(cplx(0, 1) * k * r) / (4 * pi * r);
    return {cplx(0, -1) * k * Phi, Phi * (1.0 / r - cplx(0, 1) * k) / r};
}

inline KernelValue dirac_kernel_remainder(const Vec3& z, cplx k) {
    real r = norm(z);
    if (r < 1e-12) return {k * k / (4 * pi), 0.0};
    cplx ikr = cplx(0, 1) * k * r, e = std::exp(ikr);
    cplx av = (e * (1.0 - ikr) - 1.0) / (4 * pi * r * r * r);
    cplx bs = cplx(0, -1) * k * (e - 1.0) / (4 * pi * r);
    return {bs, av};
}

}  // namespace cbem
