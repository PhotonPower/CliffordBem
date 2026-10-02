"""Eigene einfallende Felder als Python-Klassen (v0.43).

Ein einfallendes Feld ist eine Unterklasse von :class:`CustomField` mit der Methode ``fields(x)``, die fuer Punkte ``x`` der
Form (M, 3) die Felder ``(E, H)`` als komplexe Arrays der Form (M, 3) liefert -- vektorisiert mit NumPy. Normierung wie im
Kern: eps0 = mu0 = 1, exp(-i omega t), H physikalisch (ebene Welle: H = sqrt(eps/mu) d x E). Das Feld muss die Maxwell-
Gleichungen im Aussenmedium erfuellen (sonst ist die Loesung bedeutungslos); :func:`maxwell_residual` prueft das.

Der Kern ruft ``fields`` je Projektion und je Nahfeldauswertung genau einmal mit allen Punkten auf (nicht je Punkt), ausserhalb
der parallelen Schleifen; das Streufeld rechnet er ohne GIL. Verwendung::

    class Stehwelle(cb.CustomField):
        def __init__(self, medium, omega):
            self.k, self.Z = medium.k(omega), np.sqrt(medium.eps / medium.mu)
        def fields(self, x):
            c, s = np.cos(self.k * x[:, 2]), np.sin(self.k * x[:, 2])
            E = np.zeros((len(x), 3), complex); H = np.zeros_like(E)
            E[:, 0] = 2 * c; H[:, 1] = 2j * self.Z * s
            return E, H

    f = Stehwelle(water, om)
    b = f.project(P.mesh, water)                       # rechte Seite
    r = P.solve_rhs(b)
    nf = cb.exterior_near_field(P.mesh, r.h, b, water, om, f, pts)
    F = cb.force_on_sphere(cb.near_field_evaluator(P.mesh, r.h, b, water, om, f), water, (0, 0, 0), 1.5)

Felder lassen sich addieren und mit Zahlen multiplizieren, auch zusammen mit den Feldern des Kerns
(``PlaneWaveField(...) + PlaneWaveField(...)``, ``2 * f``); das ergibt ein :class:`SuperposedField`.
"""
from __future__ import annotations

from numbers import Number

import numpy as _np

from . import _cbem

__all__ = ["CustomField", "SuperposedField", "PythonPlaneWave", "AngularSpectrumField", "as_field_function", "maxwell_residual"]


def _check_pair(E, H, M):
    E = _np.asarray(E, dtype=complex)
    H = _np.asarray(H, dtype=complex)
    for name, a in (("E", E), ("H", H)):
        if a.shape != (M, 3):
            raise ValueError(f"fields(x): {name} muss die Form ({M}, 3) haben, nicht {a.shape}")
    return E, H


def _evaluate(field, x):
    """(E, H) an Punkten x (M, 3) fuer ein Feld des Kerns oder ein Python-Feld."""
    x = _np.ascontiguousarray(x, dtype=float).reshape(-1, 3)
    if isinstance(field, _cbem.IncidentField):
        E, H = field.eval(x)
        return _np.asarray(E).reshape(-1, 3), _np.asarray(H).reshape(-1, 3)
    if not hasattr(field, "fields"):
        raise TypeError("erwartet ein IncidentField des Kerns oder ein Objekt mit fields(x) -> (E, H)")
    return _check_pair(*field.fields(x), len(x))


class CustomField:
    """Basisklasse eigener einfallender Felder: ``fields(x)`` ueberschreiben.

    Attribute ``reference_E2`` und ``reference_C`` sind die Bezugsgroessen fuer Feldverstaerkung |E|^2/reference_E2 und
    optische Chiralitaet Im(E*.H)/reference_C im Nahfeld (Voreinstellung 1).
    """

    reference_E2: float = 1.0
    reference_C: float = 1.0

    def fields(self, x: _np.ndarray):
        """(E, H) an den Punkten x (M, 3); komplexe Arrays (M, 3)."""
        raise NotImplementedError("CustomField.fields(x) muss ueberschrieben werden")

    # --- Auswertung wie IncidentField des Kerns -----------------------------------------------------------------------------
    def eval(self, x):
        """(E, H) an einem Punkt (3,) oder an Punkten (M, 3)."""
        x = _np.asarray(x, dtype=float)
        E, H = _evaluate(self, x)
        return (E[0], H[0]) if x.ndim == 1 else (E, H)

    def project(self, mesh, medium, sub: int = 2):
        """Spur sqrt(eps) E + I sqrt(mu) H auf dem Netz (rechte Seite fuer solve_rhs); dieselbe Regel wie im Kern.

        Ebenes Netz (TriangleMesh): 8 Koeffizienten je Dreieck; gekruemmtes Netz (QuadraticMesh, v0.60): 24 je Element in der
        psi-Basis (wie project_incident_curved)."""
        if isinstance(mesh, _cbem.QuadraticMesh):
            pts = _cbem._quadrature_points_curved(mesh, sub)
            E, H = _evaluate(self, pts)
            return _cbem._project_samples_curved(mesh, E, H, medium, sub)
        pts = _cbem._quadrature_points(mesh, sub)
        E, H = _evaluate(self, pts)
        return _cbem._project_samples(mesh, E, H, medium, sub)

    # --- Ueberlagerung ----------------------------------------------------------------------------------------------------------
    def __add__(self, other):
        return SuperposedField([self, other])

    def __radd__(self, other):
        if isinstance(other, Number) and other == 0:              # sum([...])
            return self
        return SuperposedField([other, self])

    def __sub__(self, other):
        return SuperposedField([self, other], [1.0, -1.0])

    def __mul__(self, c):
        if not isinstance(c, Number):
            return NotImplemented
        return SuperposedField([self], [c])

    __rmul__ = __mul__

    def __neg__(self):
        return SuperposedField([self], [-1.0])

    def __repr__(self):
        return f"{type(self).__name__}()"


class SuperposedField(CustomField):
    """Lineare Ueberlagerung sum_i c_i F_i von Feldern des Kerns und Python-Feldern.

    Bezugsgroessen: wenn nicht angegeben, die der ersten Komponente (z. B. Stehwelle aus zwei ebenen Wellen: Verstaerkung
    bezogen auf eine einzelne Welle).
    """

    def __init__(self, components, coefficients=None, reference_E2=None, reference_C=None):
        comps, coefs = [], []
        coefficients = [1.0] * len(components) if coefficients is None else list(coefficients)
        if len(coefficients) != len(components):
            raise ValueError("je Komponente ein Koeffizient")
        for f, c in zip(components, coefficients):
            if isinstance(f, SuperposedField):                    # flach halten
                comps += f.components
                coefs += [c * ci for ci in f.coefficients]
            else:
                if not isinstance(f, _cbem.IncidentField) and not hasattr(f, "fields"):
                    raise TypeError(f"keine Feldkomponente: {f!r}")
                comps.append(f)
                coefs.append(complex(c))
        if not comps:
            raise ValueError("mindestens eine Komponente")
        self.components, self.coefficients = comps, coefs
        self.reference_E2 = float(comps[0].reference_E2 if reference_E2 is None else reference_E2)
        self.reference_C = float(comps[0].reference_C if reference_C is None else reference_C)

    def fields(self, x):
        E = _np.zeros((len(x), 3), complex)
        H = _np.zeros_like(E)
        for f, c in zip(self.components, self.coefficients):
            Ef, Hf = _evaluate(f, x)
            E += c * Ef
            H += c * Hf
        return E, H

    def __repr__(self):
        return f"SuperposedField({len(self.components)} Komponenten)"


class PythonPlaneWave(CustomField):
    """Ebene Welle E = p exp(i k d.x) im achiralen Medium (Referenz und Vorlage fuer eigene Felder)."""

    def __init__(self, medium, omega: float, d, p, origin=(0.0, 0.0, 0.0)):
        if medium.chiral:
            raise ValueError("PythonPlaneWave: nur achirale Medien (im chiralen Medium PlaneWaveField des Kerns verwenden)")
        self.k = medium.k(omega)
        self.d = _np.asarray(d, float) / _np.linalg.norm(d)
        self.p = _np.asarray(p, complex)
        self.origin = _np.asarray(origin, float)
        self.Z = _np.sqrt(complex(medium.eps)) / _np.sqrt(complex(medium.mu))     # wie PlaneWaveField im Kern
        self._dxp = _np.cross(self.d, self.p)
        self.reference_E2 = float(_np.sum(_np.abs(self.p) ** 2))
        self.reference_C = abs(self.Z.real) * self.reference_E2                  # |Im(E0*.H0)| der zirkularen Welle

    def fields(self, x):
        ph = _np.exp(1j * self.k * ((x - self.origin) @ self.d))[:, None]
        return ph * self.p[None, :], ph * (self.Z * self._dxp)[None, :]


class AngularSpectrumField(CustomField):
    """Ueberlagerung endlich vieler ebener Wellen (exakte Maxwell-Loesung im achiralen Medium):
    E(x) = sum_j a_j exp(i k k_j.(x - x_f)), H = sqrt(eps/mu) sum_j k_j x a_j exp(...), mit Einheitsvektoren k_j und zu k_j
    senkrechten Amplituden a_j (der Anteil parallel zu k_j wird entfernt).

    Grundlage fuer eigene Strahlen (Bessel-, Vortex-, Mehrstrahlinterferenz): siehe :meth:`bessel`.
    """

    def __init__(self, medium, omega: float, directions, amplitudes, focus=(0.0, 0.0, 0.0), reference_E2=None):
        if medium.chiral:
            raise ValueError("AngularSpectrumField: nur achirale Medien")
        kd = _np.asarray(directions, float).reshape(-1, 3)
        kd = kd / _np.linalg.norm(kd, axis=1, keepdims=True)
        a = _np.asarray(amplitudes, complex).reshape(-1, 3)
        if len(a) != len(kd):
            raise ValueError("je Richtung eine Amplitude")
        a = a - _np.sum(a * kd, axis=1, keepdims=True) * kd          # transversal
        self.k = medium.k(omega)
        self.Z = _np.sqrt(complex(medium.eps)) / _np.sqrt(complex(medium.mu))
        self.dirs, self.amps = kd, a
        self.focus = _np.asarray(focus, float)
        self.hamps = self.Z * _np.cross(kd, a)
        E0, _ = self.fields(self.focus[None, :])
        self.reference_E2 = float(_np.sum(_np.abs(E0) ** 2)) if reference_E2 is None else float(reference_E2)
        if self.reference_E2 == 0:
            self.reference_E2 = 1.0

    def fields(self, x):
        ph = _np.exp(1j * self.k * ((x - self.focus) @ self.dirs.T))        # (M, J)
        return ph @ self.amps, ph @ self.hamps

    @classmethod
    def bessel(cls, medium, omega: float, cone_angle: float, pol=(1.0, 0.0, 0.0), charge: int = 0, nphi: int = 64,
               focus=(0.0, 0.0, 0.0)):
        """Vektorieller Bessel-Strahl entlang +z: nphi ebene Wellen auf einem Kegel mit halbem Oeffnungswinkel cone_angle,
        Polarisation pol (in der Pupille, wie Richards-Wolf), topologische Ladung charge (Phase exp(i m phi): Wirbelstrahl).
        Exakte Maxwell-Loesung fuer jedes nphi; fuer nphi >> |k| rho_max gleich dem kontinuierlichen Bessel-Strahl."""
        ph = 2 * _np.pi * (_np.arange(nphi) + 0.5) / nphi
        st, ct = _np.sin(cone_angle), _np.cos(cone_angle)
        dirs = _np.column_stack([st * _np.cos(ph), st * _np.sin(ph), _np.full(nphi, ct)])
        phi_hat = _np.column_stack([-_np.sin(ph), _np.cos(ph), _np.zeros(nphi)])
        rho_hat = _np.column_stack([_np.cos(ph), _np.sin(ph), _np.zeros(nphi)])
        theta_hat = _np.column_stack([ct * _np.cos(ph), ct * _np.sin(ph), -_np.full(nphi, st)])
        pol = _np.asarray(pol, complex)
        a = (phi_hat @ pol)[:, None] * phi_hat + (rho_hat @ pol)[:, None] * theta_hat
        a *= _np.exp(1j * charge * ph)[:, None] / nphi
        return cls(medium, omega, dirs, a, focus)


def as_field_function(fun, reference_E2: float = 1.0, reference_C: float = 1.0) -> CustomField:
    """CustomField aus einer Funktion x -> (E, H)."""

    class _FunctionField(CustomField):
        def fields(self, x):
            return fun(x)

        def __repr__(self):
            return f"as_field_function({getattr(fun, '__name__', 'fun')})"

    f = _FunctionField()
    f.reference_E2, f.reference_C = float(reference_E2), float(reference_C)
    return f


def maxwell_residual(field, medium, omega: float, points, step: float | None = None) -> float:
    """Relatives Residuum der Maxwell-Gleichungen des Feldes an den Punkten (Pasteur-Medium, exp(-i omega t)):
        curl E = i omega B,  curl H = -i omega D,  D = eps E + i chi H,  B = mu H - i chi E,
    Ableitungen mit zentralen Differenzen vierter Ordnung (Schritt Voreinstellung 0,01 / |k|). Ergebnis
    ||Residuum|| / ||Terme|| ueber alle Punkte. Fuer Felder, die auf der Skala der Wellenlaenge variieren (ebene Wellen,
    Strahlen), ergibt eine Maxwell-Loesung etwa 1e-9; im Nahfeld von Quellen (Skala L < Wellenlaenge, z. B. Abstand vom
    Dipol) waechst der Differenzenfehler wie (step/L)^4 -- dann ``step`` verkleinern. Ein Feld mit falscher Normierung von H
    oder falscher Wellenzahl ergibt O(0,1..1). Funktioniert fuer Felder des Kerns und Python-Felder."""
    x = _np.ascontiguousarray(points, dtype=float).reshape(-1, 3)
    M = len(x)
    eps, mu, chi = complex(medium.eps), complex(medium.mu), complex(medium.chi)
    k = abs(omega * _np.sqrt(eps * mu)) + abs(omega * chi)
    h = 0.01 / k if step is None else float(step)
    offs = [(-2, 1.0), (-1, -8.0), (1, 8.0), (2, -1.0)]                      # (f(x-2h) - 8 f(x-h) + 8 f(x+h) - f(x+2h)) / 12h
    stencil = [x]
    for i in range(3):
        for o, _ in offs:
            e = _np.zeros(3)
            e[i] = o * h
            stencil.append(x + e)
    E, H = _evaluate(field, _np.vstack(stencil))
    E0, H0 = E[:M], H[:M]
    dE = _np.zeros((M, 3, 3), complex)                                      # dE[:, i, j] = d_i E_j
    dH = _np.zeros_like(dE)
    blk = 1
    for i in range(3):
        for o, w in offs:
            sl = slice(blk * M, (blk + 1) * M)
            dE[:, i, :] += w * E[sl] / (12 * h)
            dH[:, i, :] += w * H[sl] / (12 * h)
            blk += 1

    def curl(d):
        return _np.stack([d[:, 1, 2] - d[:, 2, 1], d[:, 2, 0] - d[:, 0, 2], d[:, 0, 1] - d[:, 1, 0]], axis=1)

    D = eps * E0 + 1j * chi * H0
    B = mu * H0 - 1j * chi * E0
    rE = curl(dE) - 1j * omega * B
    rH = curl(dH) + 1j * omega * D
    scale = _np.linalg.norm(omega * B) + _np.linalg.norm(omega * D)
    if scale == 0:
        raise ValueError("maxwell_residual: Feld ist an allen Punkten null")
    return float((_np.linalg.norm(rE) + _np.linalg.norm(rH)) / scale)
