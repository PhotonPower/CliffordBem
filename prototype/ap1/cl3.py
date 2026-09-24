"""
Minimale Implementierung der komplexifizierten Clifford-Algebra Cl_3(C).

Basisblades werden als Bitmasken kodiert: e1=0b001, e2=0b010, e3=0b100,
e12=0b011, e13=0b101, e23=0b110, e123=0b111. Ein Multivektor ist ein
Array der Laenge 8 (numpy, dtype complex oder object fuer SymPy).
Die imaginaere Einheit i ist der Skalar 1j bzw. sympy.I und kommutiert mit
allem; der Pseudoskalar I = e123 ist davon strikt verschieden.
"""
import numpy as np

GRADE = np.array([bin(b).count("1") for b in range(8)])
NAMES = ["1", "e1", "e2", "e12", "e3", "e13", "e23", "e123"]


def _reorder_sign(a, b):
    a >>= 1
    s = 0
    while a:
        s += bin(a & b).count("1")
        a >>= 1
    return -1 if (s & 1) else 1


# Multiplikationstabelle: blade a * blade b = SIGN[a,b] * blade (a ^ b)
SIGN = np.array([[_reorder_sign(a, b) for b in range(8)] for a in range(8)])


class MV:
    """Multivektor in Cl_3(C). Koeffizienten numerisch oder symbolisch."""

    def __init__(self, c=None, symbolic=False):
        if c is None:
            c = [0] * 8
        self.c = np.array(c, dtype=object if symbolic else complex)
        self.symbolic = symbolic

    # --- Konstruktion ---
    @staticmethod
    def blade(b, coef=1, symbolic=False):
        c = [0] * 8
        c[b] = coef
        return MV(c, symbolic)

    @staticmethod
    def vec(v, symbolic=False):
        c = [0] * 8
        c[1], c[2], c[4] = v[0], v[1], v[2]
        return MV(c, symbolic)

    # --- Arithmetik ---
    def _coerce(self, o):
        if isinstance(o, MV):
            return o
        return MV.blade(0, o, self.symbolic)

    def __add__(self, o):
        o = self._coerce(o)
        return MV(self.c + o.c, self.symbolic or o.symbolic)

    __radd__ = __add__

    def __neg__(self):
        return MV(-self.c, self.symbolic)

    def __sub__(self, o):
        return self + (-self._coerce(o))

    def __rsub__(self, o):
        return self._coerce(o) - self

    def __mul__(self, o):
        if not isinstance(o, MV):
            return MV(self.c * o, self.symbolic)
        sym = self.symbolic or o.symbolic
        r = [0] * 8
        for a in range(8):
            if self.c[a] == 0:
                continue
            for b in range(8):
                if o.c[b] == 0:
                    continue
                r[a ^ b] = r[a ^ b] + SIGN[a, b] * self.c[a] * o.c[b]
        return MV(r, sym)

    def __rmul__(self, o):  # Skalar * MV
        return MV(self.c * o, self.symbolic)

    # --- Involutionen und Projektionen ---
    def grade(self, *ks):
        c = [self.c[b] if GRADE[b] in ks else 0 for b in range(8)]
        return MV(c, self.symbolic)

    def reverse(self):
        s = [1, 1, -1, -1]  # (-1)^{k(k-1)/2}
        return MV([s[GRADE[b]] * self.c[b] for b in range(8)], self.symbolic)

    def involute(self):
        return MV([(-1) ** GRADE[b] * self.c[b] for b in range(8)], self.symbolic)

    def conj(self):
        if self.symbolic:
            import sympy as sp
            return MV([sp.conjugate(x) for x in self.c], True)
        return MV(np.conj(self.c))

    def scalar(self):
        return self.c[0]

    def norm2(self):
        """<A ~A^*>_0 = Summe |a_B|^2 (positiv definit)."""
        return (self * self.reverse().conj()).scalar()

    def simplify(self):
        import sympy as sp
        return MV([sp.simplify(x) for x in self.c], True)

    def is_zero(self, tol=1e-10):
        if self.symbolic:
            import sympy as sp
            return all(sp.simplify(x) == 0 for x in self.c)
        return np.max(np.abs(self.c)) < tol

    def __repr__(self):
        parts = [f"({self.c[b]})*{NAMES[b]}" for b in range(8) if self.c[b] != 0]
        return " + ".join(parts) if parts else "0"

    # --- Matrixdarstellung der Links-/Rechtsmultiplikation ---
    def left_matrix(self):
        M = np.zeros((8, 8), dtype=complex)
        for b in range(8):
            M[:, b] = (self * MV.blade(b)).c
        return M

    def right_matrix(self):
        M = np.zeros((8, 8), dtype=complex)
        for b in range(8):
            M[:, b] = (MV.blade(b) * self).c
        return M


E1, E2, E3 = MV.blade(1), MV.blade(2), MV.blade(4)
PSEUDO = MV.blade(7)  # I = e1 e2 e3


# --- Vektorableitung (symbolisch) ---
def nabla_left(F, X):
    """nabla F = sum_j e_j d_j F   (X = (x, y, z) SymPy-Symbole)."""
    import sympy as sp
    r = MV(symbolic=True)
    for j, b in enumerate([1, 2, 4]):
        dF = MV([sp.diff(x, X[j]) for x in F.c], True)
        r = r + MV.blade(b, 1, True) * dF
    return r


def nabla_right(F, X):
    """F nabla(von rechts) = sum_j (d_j F) e_j."""
    import sympy as sp
    r = MV(symbolic=True)
    for j, b in enumerate([1, 2, 4]):
        dF = MV([sp.diff(x, X[j]) for x in F.c], True)
        r = r + dF * MV.blade(b, 1, True)
    return r
