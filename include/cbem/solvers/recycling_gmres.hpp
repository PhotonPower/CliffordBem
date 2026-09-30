#pragma once
// GMRES mit Krylov-Recycling (GCRO-Prinzip) fuer Folgen von Gleichungssystemen mit derselben Matrix und verschiedenen rechten
// Seiten (v0.28), etwa die Orientierungen eines Dipols oder viele Emitterorte auf einem Netz.
//
// Rechtsvorkonditioniert wie gmres: B = A M, geloest wird B z = b, x = M z. Gespeichert werden U (im z-Raum) und C = B U mit
// orthonormalem C. Fuer eine neue rechte Seite:
//   1. Projektion: z0 = U C^H b, r0 = b - C C^H b;
//   2. GMRES auf dem deflationierten Operator (I - C C^H) B: jeder neue Arnoldi-Vektor wird zuerst gegen C orthogonalisiert
//      (zweifaches klassisches Gram-Schmidt), B V_m = C E + V_{m+1} Hq;
//   3. y = argmin |beta e1 - Hq y|, z = z0 + (V_m - U E) y (exakt, da B (V_m - U E) y = V_{m+1} Hq y);
//   4. Erweiterung des Unterraums: QR von Hq = Q R, C_neu = V_{m+1} Q, U_neu = (V_m - U E) R^{-1}, angehaengt bis max_recycle.
// Die Richtungen, die fruehere Loesungen gefunden haben, muessen nicht neu gesucht werden. Neustarts innerhalb einer Loesung
// erweitern den Unterraum ebenfalls. Speicher: 2 max_recycle Vektoren der Laenge n.
#include <vector>
#include "cbem/solvers/gmres.hpp"

namespace cbem {

class RecyclingGmres {
public:
    RecyclingGmres(LinOp A, const LinOp* M, std::size_t max_recycle = 120) : A_(std::move(A)), max_recycle_(max_recycle) {
        if (M) { M_ = *M; hasM_ = true; }
    }
    // Loest A x = b (Start bei x = 0); tol relativ zu |b|
    GmresResult solve(const std::vector<cplx>& b, std::vector<cplx>& x, real tol = 1e-6, int restart = 200, int max_iter = 2000);
    std::size_t recycled() const { return U_.size(); }
    void clear() { U_.clear(); C_.clear(); }
private:
    void applyB(const std::vector<cplx>& v, std::vector<cplx>& out);
    LinOp A_, M_;
    bool hasM_ = false;
    std::size_t max_recycle_;
    std::vector<std::vector<cplx>> U_, C_;
    std::vector<cplx> tmp_;
};

}  // namespace cbem
