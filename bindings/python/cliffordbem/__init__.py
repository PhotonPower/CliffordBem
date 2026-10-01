"""CliffordBem -- koordinatenfreie Randelementmethode der Nano-Optik in Cl3(C).

Python-Anbindung (v0.42) an den C++-Kern. Alle Klassen und Funktionen des Erweiterungsmoduls ``_cbem`` sind
hier direkt verfuegbar; dazu kommen einige Komfortfunktionen (Materialien, Wellenlaengen, Spektren).

Einheiten wie im Kern: Laengen in einer frei gewaehlten Einheit L (z. B. ``unit = 20`` nm), eps0 = mu0 = 1,
omega = k0 = 2 pi L / lambda, Zeitkonvention exp(-i omega t). Spuren sind flache complex128-Arrays der Laenge
8 N (Blade-Reihenfolge 1, e1, e2, e12, e3, e13, e23, e123 je Dreieck).

Beispiel (Goldkugel mit 40 nm Durchmesser in Wasser bei 530 nm)::

    import cliffordbem as cb
    gold = cb.make_material("Au")
    unit, lam = 20.0, 530.0
    P = cb.ScatteringProblem(cb.make_icosphere(8), gold.medium(lam), cb.omega_from_wavelength(lam, unit),
                             outer=cb.Medium(eps=1.33**2))
    r = P.solve_plane_wave(d=(0, 0, 1), p=(1, 0, 0))
    print(r.sigma_ext * unit**2, "nm^2")
"""
from __future__ import annotations

import os as _os

import numpy as _np

from . import _cbem
from ._cbem import *  # noqa: F401,F403
from ._cbem import __version__, _make_material  # noqa: F401  __version__ wird re-exportiert

__all__ = [n for n in dir(_cbem) if not n.startswith("_")] + [
    "__version__", "DATA_DIR", "make_material", "omega_from_wavelength", "wavelength_from_omega", "trace_blocks", "polarization_basis",
    "spectrum", "e1", "e2", "e3", "e12", "e13", "e23", "I",
]


def _find_data_dir() -> str:
    """Materialdaten: CBEM_DATA_DIR, sonst im Paket (Installation), sonst im Quellbaum (Bau mit CMake)."""
    env = _os.environ.get("CBEM_DATA_DIR")
    if env:
        return env
    here = _os.path.dirname(_os.path.abspath(__file__))
    pkg = _os.path.join(here, "data", "materials")
    if _os.path.isdir(pkg):
        return pkg
    src = getattr(_cbem, "_SOURCE_DATA_DIR", "")
    if src and _os.path.isdir(src):
        return src
    return _os.path.join("data", "materials")


#: Verzeichnis der Materialtabellen (Johnson-Christy Au, Ag)
DATA_DIR: str = _find_data_dir()

# Basis-Blades (unveraenderliche Multivektoren)
e1 = Multivector.blade("e1")      # noqa: F405
e2 = Multivector.blade("e2")      # noqa: F405
e3 = Multivector.blade("e3")      # noqa: F405
e12 = Multivector.blade("e12")    # noqa: F405
e13 = Multivector.blade("e13")    # noqa: F405
e23 = Multivector.blade("e23")    # noqa: F405
I = Multivector.blade("e123")     # noqa: F405,E741  Pseudoskalar


def make_material(spec: str, data_dir: str | None = None):
    """Material aus einer Kurzangabe: ``"Au"``/``"Ag"`` (Johnson-Christy), ein Pfad zu einer
    refractiveindex.info-YAML-Datei oder eine konstante Permittivitaet ``"re,im"``."""
    if spec.endswith((".yml", ".yaml")) and _os.path.isfile(spec):
        return TabulatedMaterial.from_yaml(spec)  # noqa: F405
    return _make_material(spec, data_dir or DATA_DIR)


def omega_from_wavelength(lambda_nm, unit_nm: float):
    """omega = k0 = 2 pi unit / lambda (Laengen in Einheiten von unit_nm); Zahl oder Array."""
    w = 2.0 * _np.pi * unit_nm / _np.asarray(lambda_nm, dtype=float)
    return float(w) if w.ndim == 0 else w


def wavelength_from_omega(omega, unit_nm: float):
    """Umkehrung von omega_from_wavelength: lambda in nm."""
    lam = 2.0 * _np.pi * unit_nm / _np.asarray(omega, dtype=float)
    return float(lam) if lam.ndim == 0 else lam


def trace_blocks(h) -> _np.ndarray:
    """Spur (8 N,) als (N, 8): Multivektor-Koeffizienten je Dreieck."""
    h = _np.asarray(h)
    if h.size % 8:
        raise ValueError("Laenge der Spur ist kein Vielfaches von 8")
    return h.reshape(-1, 8)


def polarization_basis(d):
    """Orthonormale reelle Polarisationen (u, v) zur Richtung d, (u, v, d) rechtshaendig, wie circular_polarization."""
    c = circular_polarization(d, +1)  # noqa: F405  c = u + i v
    return _np.real(c), _np.imag(c)


def spectrum(bodies, materials, wavelengths, unit: float, n_bg: float = 1.0, pol: str = "lin", orient: int = 1,
             chi=None, host_chi: complex = 0.0, hmatrix=None, options=None, callback=None) -> dict:
    """Extinktionsspektrum homogener Koerper wie die App ``spectrum`` (ohne Beschichtungen).

    bodies: TriangleMesh oder Liste davon (Laengen in Einheiten von ``unit`` nm); materials: Material oder Liste je
    Koerper (z. B. ``make_material("Au")``); wavelengths in nm; ``pol``: ``"lin"`` (bei orient > 1 unpolarisiert
    gemittelt) oder ``"circ"`` (sigma+, sigma-, CD); ``orient``: Lebedev-Punkte 1, 6, 14, 26 (Einfallsrichtungen);
    ``chi``: Chiralitaet je Koerper; ``host_chi``: chirales Aussenmedium (nur ``pol="circ"``).
    ``callback(i, row)`` wird nach jeder Wellenlaenge aufgerufen (Fortschritt).

    Gibt ein dict von Arrays zurueck: lambda_nm, sigma_nm2, sigma_plus_nm2, sigma_minus_nm2, cd_nm2,
    forward (S(0), komplex), iterations.
    """
    if isinstance(bodies, TriangleMesh):  # noqa: F405
        bodies = [bodies]
    bodies = list(bodies)
    if not isinstance(materials, (list, tuple)):
        materials = [materials] * len(bodies)
    if len(materials) != len(bodies):
        raise ValueError("je Koerper ein Material erwartet")
    chi = [0.0] * len(bodies) if chi is None else (list(chi) if _np.ndim(chi) else [chi] * len(bodies))
    if pol not in ("lin", "circ"):
        raise ValueError("pol: 'lin' oder 'circ'")
    if host_chi != 0 and pol != "circ":
        raise ValueError("chirales Aussenmedium nur mit zirkularer Polarisation (pol='circ')")
    dirs, weights = lebedev(orient)  # noqa: F405
    hp = hmatrix if hmatrix is not None else HMatrixParams()  # noqa: F405
    so = options if options is not None else SolveOptions()  # noqa: F405
    bg = Medium(eps=n_bg * n_bg, chi=host_chi)  # noqa: F405
    lams = _np.atleast_1d(_np.asarray(wavelengths, dtype=float))
    out = {k: _np.zeros(len(lams)) for k in ("sigma_nm2", "sigma_plus_nm2", "sigma_minus_nm2", "cd_nm2")}
    out["forward"] = _np.zeros(len(lams), dtype=complex)
    out["iterations"] = _np.zeros(len(lams), dtype=int)
    u2 = unit * unit
    for i, lam in enumerate(lams):
        om = omega_from_wavelength(lam, unit)
        media = [Medium(eps=m.eps(lam), chi=c) for m, c in zip(materials, chi)]  # noqa: F405
        P = ScatteringProblem(bodies, media, om, outer=bg, hmatrix=hp)  # noqa: F405
        s = sp = sm = 0.0
        S = 0j
        its = 0
        for d, w in zip(dirs, weights):
            if pol == "circ":
                a = P.solve_plane_wave(d, circular_polarization(d, +1), so)  # noqa: F405
                b = P.solve_plane_wave(d, circular_polarization(d, -1), so)  # noqa: F405
                sp += w * a.sigma_ext
                sm += w * b.sigma_ext
                S += w * 0.5 * (a.forward + b.forward)
                its = max(its, a.iterations, b.iterations)
            else:
                u, v = polarization_basis(d)
                a = P.solve_plane_wave(d, u, so)
                sa, Sa = a.sigma_ext, a.forward
                its = max(its, a.iterations)
                if len(dirs) > 1:
                    b = P.solve_plane_wave(d, v, so)
                    sa, Sa = 0.5 * (sa + b.sigma_ext), 0.5 * (Sa + b.forward)
                    its = max(its, b.iterations)
                s += w * sa
                S += w * Sa
        if pol == "circ":
            s = 0.5 * (sp + sm)
        out["sigma_nm2"][i] = s * u2
        out["sigma_plus_nm2"][i] = sp * u2
        out["sigma_minus_nm2"][i] = sm * u2
        out["cd_nm2"][i] = (sp - sm) * u2
        out["forward"][i] = S
        out["iterations"][i] = its
        if callback is not None:
            callback(i, {k: v[i] for k, v in out.items()} | {"lambda_nm": lam})
    out["lambda_nm"] = lams
    return out
