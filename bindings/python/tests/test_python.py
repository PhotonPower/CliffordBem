"""Tests der Python-Anbindung (v0.42/v0.43): Algebra, Netze, Materialien, Parameter, Streuung (Regression gegen den C++-Test und
bitgleich zum Low-Level-Pfad), Mie, chirale Medien, mehrere Koerper, Nahfeld gegen Mie, Kraefte (drei Wege und Mie), Dipol,
Strahlen, geschichtete Probleme, Fehlerbehandlung, Lebensdauer (keep_alive), Threads ohne GIL; eigene einfallende Felder
(bitgleich zum Kern, Ueberlagerung, Maxwell-Pruefung, Bessel-Strahlen, Fehlerbehandlung).

Aufruf: PYTHONPATH=build/python python3 bindings/python/tests/test_python.py   (oder pytest; ctest: test_python)
Referenzen aus tools/ (mie.py, mie_nearfield.py, mie_force.py): Verzeichnis ueber CBEM_TOOLS_DIR, sonst relativ zum Quellbaum.
"""
import copy
import os
import pickle
import sys
import tempfile
import threading
import time
import traceback

import numpy as np

import cliffordbem as cb

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("CBEM_TOOLS_DIR", os.path.join(_here, "..", "..", "..", "tools")))
import mie  # noqa: E402
import mie_force  # noqa: E402
import mie_nearfield  # noqa: E402

Z, X = (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)
GOLD = cb.Medium(eps=-11 + 1.2j)
WATER = cb.Medium(eps=1.7689)


def raises(exc, f, *a, **k):
    try:
        f(*a, **k)
    except exc:
        return True
    return False


# --- Algebra --------------------------------------------------------------------------------------------------------------------
def test_algebra():
    e1, e2, e3, I = cb.e1, cb.e2, cb.e3, cb.I
    assert (e1 * e2)["e12"] == 1 and (e2 * e1)["e12"] == -1
    assert (I * I) == cb.Multivector.scalar(-1)
    for e in (e1, e2, e3, cb.e12):                                     # I ist zentral
        assert I * e == e * I
    A = cb.Multivector(np.arange(8) + 1j * np.arange(8)[::-1])
    B = cb.Multivector(np.cos(np.arange(8)) + 0.5j)
    assert np.allclose((A * B).c, A.left_matrix() @ B.c)                 # Linksmultiplikation
    assert np.allclose((A * A.inverse()).c, cb.Multivector.scalar(1).c)
    P = (1 + 1j * I) / 2                                                 # Helizitaetsprojektor: idempotent, Nullteiler
    assert np.allclose((P * P).c, P.c)
    assert np.allclose(P.left_matrix(), cb.helicity_projector(+1))
    assert raises(ValueError, P.inverse)
    assert (A.reverse().reverse()) == A and A.grade(2)["e12"] == A["e12"] and A.grade(2)["e1"] == 0
    assert pickle.loads(pickle.dumps(A)) == A and copy.deepcopy(A) == A
    assert "e12" in repr(e1 * e2) and raises(KeyError, lambda: A["e4"])


# --- Netze ------------------------------------------------------------------------------------------------------------------------
def test_mesh():
    m = cb.make_icosphere(4)
    assert len(m) == 320 and m.unknowns == 2560 and m.points.shape == (162, 3) and m.triangles.shape == (320, 3)
    assert np.allclose(np.linalg.norm(m.normals, axis=1), 1) and abs(m.areas.sum() - 4 * np.pi) < 0.3
    assert np.all(np.einsum("ij,ij->i", m.normals, m.centroids) > 0)    # nach aussen
    m2 = cb.TriangleMesh(m.points, m.triangles)
    assert np.allclose(m2.areas, m.areas) and np.allclose(pickle.loads(pickle.dumps(m)).normals, m.normals)
    flipped = cb.TriangleMesh(m.points, m.triangles[:, ::-1], orient_outward=True)
    assert np.allclose(flipped.normals, m.normals)
    t = cb.translated(m, (3, 0, 0), 2.0)
    assert np.allclose(t.points.mean(axis=0), [3, 0, 0], atol=1e-12) and abs(cb.signed_volume(t) / cb.signed_volume(m) - 8) < 1e-12
    o = cb.offset_surface(m, 0.1)
    assert abs(np.linalg.norm(o.points, axis=1).mean() - 1.1) < 0.01
    w = cb.winding_number(m, np.array([[0, 0, 0], [2, 0, 0]]))
    assert abs(w[0] - 1) < 1e-10 and abs(w[1]) < 1e-10 and abs(cb.winding_number(m, (0.1, 0.2, 0)) - 1) < 1e-10
    d = cb.distance_to_surface(m, np.array([[0, 0, 2.0]]), 5.0)
    assert abs(d[0] - 1) < 0.02
    assert raises(IndexError, cb.TriangleMesh, m.points, m.triangles + 1000)
    assert raises(ValueError, cb.TriangleMesh, m.points[:, :2], m.triangles)
    assert raises(RuntimeError, cb.require_separated, m, cb.translated(m, (1.0, 0, 0)), 0.0, "Test")
    with tempfile.TemporaryDirectory() as tmp:                           # Gmsh-Rundreise mit zwei Koerpern
        f = os.path.join(tmp, "two.msh")
        cb.write_gmsh22(f, [m, cb.translated(m, (3, 0, 0))])
        bodies = cb.read_gmsh(f)
        assert [tag for tag, _ in bodies] == [1, 2] and np.allclose(bodies[1][1].areas.sum(), m.areas.sum())
    dirs, wts = cb.lebedev(14)
    assert dirs.shape == (14, 3) and abs(wts.sum() - 1) < 1e-12


# --- Materialien, Medien, Parameter ---------------------------------------------------------------------------------------------
def test_materials_and_params():
    au = cb.make_material("Au")
    e = au.eps(530.0)
    assert -5.5 < e.real < -4 and 1.5 < e.imag < 3                      # Johnson-Christy um 530 nm
    assert np.allclose(au.eps(np.array([500.0, 530.0]))[1], e) and au.medium(530.0).eps == e
    assert cb.make_material("2.25,0.1").eps(600) == 2.25 + 0.1j and cb.ConstantMaterial(4).eps(1) == 4
    tab = cb.TabulatedMaterial("x", [400, 800], [1.5, 1.6], [0, 0.2])
    assert np.isclose(tab.eps(600), (1.55 + 0.1j) ** 2)
    m = cb.Medium(eps=2.25, chi=0.1)
    assert np.isclose(m.k(1.0, +1), 1.4) and np.isclose(m.k(1.0, -1), 1.6) and m.chiral and not GOLD.chiral
    hp = cb.HMatrixParams(eps=1e-5, leaf=16)
    assert hp.eps == 1e-5 and hp.leaf == 16 and hp.eta == 1.0 and "leaf=16" in repr(hp)
    assert raises(AttributeError, cb.HMatrixParams, epsilon=1)
    assert cb.SolveOptions(tol=1e-9).copy().tol == 1e-9
    assert np.isclose(cb.omega_from_wavelength(2 * np.pi * 20, 20), 1.0)
    assert np.allclose(cb.wavelength_from_omega(cb.omega_from_wavelength([400, 500], 20), 20), [400, 500])


# --- Streuung --------------------------------------------------------------------------------------------------------------------
def _glass_regression_parts():
    m = cb.make_icosphere(4)
    ep = cb.EntryParams(sauter_schwab=False, adaptive_outer=False, near_subdivision=4)   # Nahfeldregel des Python-Prototyps
    return m, ep, cb.HMatrixParams(eps=1e-8), cb.SolveOptions(tol=1e-10, restart=300, max_iter=1000)


def test_scattering_regression_and_low_level():
    m, ep, hp, so = _glass_regression_parts()
    P = cb.ScatteringProblem(m, cb.Medium(eps=2.25), 1.0, hmatrix=hp, entries=ep)
    r = P.solve_plane_wave(Z, X, so)
    assert abs(r.sigma_ext / np.pi - 0.20291) < 2e-5, r                  # wie tests/test_scattering.cpp
    # derselbe Weg zu Fuss (wie apps/scatter_sphere.cpp); Zwischenobjekte werden freigegeben, keep_alive haelt die Kette
    Ein, Eout = cb.KernelEntries(m, 1.5, ep), cb.KernelEntries(m, 1.0, ep)
    E1, E2 = cb.CauchyOperator(m, cb.KernelHMatrix(Ein, hp)), cb.CauchyOperator(m, cb.KernelHMatrix(Eout, hp))
    del Ein, Eout
    T = cb.TransmissionOperator(m, E1, E2, cb.Medium(eps=2.25), cb.Medium())
    del E1, E2, m
    import gc
    gc.collect()
    mesh = P.mesh
    b = cb.project_plane_wave(mesh, 1.0, 1.0, Z, X)
    h, g = cb.gmres(T, b, M=T.precondition, tol=1e-10, restart=300, max_iter=1000)
    assert g.converged and g.iterations == r.iterations and np.array_equal(h, r.h)
    assert np.allclose(T.apply(h), P.apply(h)) and np.allclose(T.apply(h), b, atol=1e-8)
    # Plemelj: E ist eine Involution (E E = 1) bis auf Diskretisierung/Kompression
    hs = h - b
    assert np.isclose(cb.extinction_cross_section(mesh, hs, 1.0, 1.0, Z, X), r.sigma_ext)
    assert np.isclose(cb.forward_amplitude(mesh, hs, 1.0, 1.0, Z, X), r.forward)
    assert np.isclose(r.sigma_ext, 4 * np.pi * r.forward.real)           # sigma_ext = 4 pi / k^2 Re S(0)
    Einf = cb.far_field_E(mesh, hs, 1.0, 1.0, np.array([Z]))[0]
    assert np.isclose(4 * np.pi * np.imag(Einf[0]), r.sigma_ext)          # optisches Theorem
    # beliebige rechte Seite = ebene Welle
    r2 = P.solve_rhs(cb.plane_wave_trace(mesh, cb.Medium(), 1.0, Z, X), so)
    assert np.allclose(r2.h, r.h)


def test_gold_against_mie():
    q_mie = mie.qext(0.5, -11 + 1.2j)
    P = cb.ScatteringProblem(cb.make_icosphere(4), GOLD, 0.5, hmatrix=cb.HMatrixParams(eps=1e-8))
    q = P.solve_plane_wave(Z, X, cb.SolveOptions(tol=1e-10)).sigma_ext / np.pi
    assert abs(q - q_mie) < abs(0.66915 - q_mie), (q, q_mie)             # Sauter-Schwab naeher an Mie (wie test_scattering)
    # Hintergrundmedium: Kugel aus dem Aussenmedium ist unsichtbar
    P0 = cb.ScatteringProblem(cb.make_icosphere(4), WATER, 0.5, outer=WATER)
    assert abs(P0.solve_plane_wave(Z, X).sigma_ext) < 1e-6


def test_chiral():
    m = cb.make_icosphere(4)
    host = cb.Medium(eps=1.7689, chi=0.02)
    P = cb.ScatteringProblem(m, host, 0.5, outer=host)                     # unsichtbar auch im chiralen Medium
    assert abs(P.solve_plane_wave(Z, cb.circular_polarization(Z, +1)).sigma_ext) < 1e-6
    assert raises(ValueError, P.solve_plane_wave, Z, X)                   # linear im chiralen Medium abgelehnt
    # chirale Kugel: Spiegelsymmetrie sigma_s(chi) = sigma_-s(-chi)
    so = cb.SolveOptions(tol=1e-9)
    s = {}
    for chi in (0.1, -0.1):
        Q = cb.ScatteringProblem(m, cb.Medium(eps=2.25, chi=chi), 1.0)
        for sg in (+1, -1):
            s[chi, sg] = Q.solve_plane_wave(Z, cb.circular_polarization(Z, sg), so).sigma_ext
    assert abs(s[0.1, +1] - s[-0.1, -1]) < 1e-5 * s[0.1, 1] and abs(s[0.1, 1] - s[0.1, -1]) > 1e-4 * s[0.1, 1]


def test_multibody_and_spectrum():
    m = cb.make_icosphere(3)
    single = cb.ScatteringProblem(m, GOLD, 0.5).solve_plane_wave(Z, X).sigma_ext
    P = cb.ScatteringProblem([m, cb.translated(m, (30, 0, 0))], [GOLD, GOLD], 0.5)
    assert P.body_begin == [0, 180, 360] and len(P.bodies) == 2
    assert abs(P.solve_plane_wave(Z, X).sigma_ext / (2 * single) - 1) < 0.05       # weit getrennt: additiv
    assert raises(ValueError, cb.ScatteringProblem, [m, m], [GOLD], 0.5)
    # spectrum() wie die App (eine Wellenlaenge, direkte Rechnung zum Vergleich)
    au = cb.make_material("Au")
    sp = cb.spectrum(m, au, [520.0], unit=20.0, n_bg=1.33)
    om = cb.omega_from_wavelength(520.0, 20.0)
    ref = cb.ScatteringProblem(m, au.medium(520.0), om, outer=cb.Medium(eps=1.33 ** 2)).solve_plane_wave(Z, X)
    assert np.isclose(sp["sigma_nm2"][0], ref.sigma_ext * 400) and sp["lambda_nm"][0] == 520.0
    spc = cb.spectrum(m, au, [520.0], unit=20.0, n_bg=1.33, pol="circ")
    assert abs(spc["cd_nm2"][0]) < 1e-6 * spc["sigma_nm2"][0]            # achirale Kugel: kein CD


# --- Nahfeld und Kraefte ----------------------------------------------------------------------------------------------------------
def test_near_field_against_mie():
    om, n = 0.5, 6
    P = cb.ScatteringProblem(cb.make_icosphere(n), GOLD, om, outer=WATER)
    r = P.solve_plane_wave(Z, X, cb.SolveOptions(tol=1e-8))
    pts = np.array([[1.5, 0, 0], [0, 1.4, 0.3], [0.2, -0.3, -1.6], [2.5, 1.0, 1.0], [0.9, 0.9, 0]])
    nf = cb.exterior_near_field(P.mesh, r.h, WATER, om, Z, X, pts)
    assert nf["E"].shape == (5, 3) and not nf["inside"].any()
    for i, x in enumerate(pts):
        E, H = mie_nearfield.fields(om, [1.0], [-11 + 1.2j], 1.7689, x)
        # Fehler bezogen auf max(|E|, |E0|): in Feldminima ist der relative Fehler kein sinnvolles Mass; er faellt mit O(h^2)
        # (Punkt 1: 0,040 / 0,018 / 0,010 fuer n = 4 / 6 / 8)
        assert np.linalg.norm(nf["E"][i] - E) < 0.03 * max(np.linalg.norm(E), 1.0), (i, nf["E"][i], E)
        assert np.linalg.norm(nf["H"][i] - H) < 0.03 * max(np.linalg.norm(H), np.sqrt(1.7689))
    assert np.allclose(nf["enhancement"], np.sum(np.abs(nf["E"]) ** 2, axis=1))       # |E0| = 1
    inside = cb.exterior_near_field(P.mesh, r.h, WATER, om, Z, X, np.array([0.0, 0.0, 0.0]))
    assert inside["inside"][0]
    # Auswerter = exterior_near_field, haelt eigene Kopien (Spur und Problem duerfen verschwinden)
    ev = cb.plane_wave_evaluator(P.mesh, r.h, WATER, om, Z, X)
    h_copy = r.h
    del P, r
    nf2 = ev(pts)
    assert np.allclose(nf2["E"], nf["E"])
    # Nahfeldoperator (H-Matrix) = direkte Summation
    m = cb.make_icosphere(n)
    b = cb.plane_wave_trace(m, WATER, om, Z, X)
    k = WATER.k(om)
    far_pts = np.random.default_rng(1).normal(size=(50, 3))
    far_pts = 3 * far_pts / np.linalg.norm(far_pts, axis=1, keepdims=True)
    direct = cb.scattered_field(m, h_copy - b, k, far_pts)
    op = cb.NearFieldOperator(m, k, far_pts, cb.HMatrixParams(eps=1e-8))
    del m
    assert np.allclose(op.apply(h_copy - b), direct, atol=1e-6 * np.abs(direct).max())


def test_forces():
    om, n = 0.5, 6
    P = cb.ScatteringProblem(cb.make_icosphere(n), GOLD, om, outer=WATER)
    r = P.solve_plane_wave(Z, X, cb.SolveOptions(tol=1e-9))
    Ft = cb.force_from_traces(P.mesh, r.h, WATER, P.body_begin)[0]
    Fs = cb.force_on_sphere(P.mesh, r.h, WATER, om, Z, X, (0, 0, 0), 1.5)
    Ff = cb.force_from_far_field(P.mesh, r.h, WATER, om, Z, X, r.sigma_ext)
    ev = cb.plane_wave_evaluator(P.mesh, r.h, WATER, om, Z, X)
    Fe = cb.force_on_sphere(ev, WATER, (0, 0, 0), 1.5)
    assert np.allclose(Fe, Fs)
    for F in (Fs, Ff):
        assert abs(F[2] / Ft[2] - 1) < 0.01 and abs(F[0]) < 1e-3 * abs(F[2])   # drei Wege gleich (test_optical_force: 0,3 %)
    Fmie = mie_force.force_z(om, [1.0], [-11 + 1.2j], 1.7689)
    assert abs(Fs[2] / Fmie - 1) < 0.05, (Fs, Fmie)
    # Dipolteilchen in der ebenen Welle: Strahlungsdruck = Mie mit n = 1 (nur Streuung am Teilchen, Koerper weit weg ist egal)
    g = cb.fields_with_gradients(ev, np.array([[4.0, 0, 0]]), 1e-3)[0]
    assert g.dE.shape == (3, 3) and np.all(np.isfinite(g.E))
    a = cb.DipolePolarizability(Ae=1e-3 + 1e-3j)
    assert np.all(np.isfinite(cb.dipole_particle_force(g, a, WATER, om)))


def test_dipole_and_beams():
    om = 0.5
    # Dipol ohne Streuer (Kugel aus dem Aussenmedium): Raten = 1 (wie test_dipole)
    m = cb.make_icosphere(4)
    P = cb.ScatteringProblem(m, WATER, om, outer=WATER)
    x0, p = (0, 0, 1.3), (0, 0, 1)
    b = cb.project_dipole(m, WATER, om, x0, p)
    r = P.solve_rhs(b, cb.SolveOptions(tol=1e-10))
    R = cb.dipole_rates(P.mesh, r.h, b, WATER, om, x0, p)
    assert abs(R.total - 1) < 1e-6 and abs(R.radiative - 1) < 1e-5 and abs(R.quantum_yield - 1) < 1e-5
    E, H = cb.dipole_field(WATER, om, x0, p, np.array([[0, 0, 5.0], [3, 0, 1.3]]))
    E1, H1 = cb.DipoleField(WATER, om, x0, p).eval(np.array([3, 0, 1.3]))
    assert E.shape == (2, 3) and np.allclose(E[1], E1) and np.allclose(H[1], H1)
    assert np.allclose(cb.DipoleField(WATER, om, x0, p).project(m, WATER), cb.IncidentField.project(cb.DipoleField(WATER, om, x0, p), m, WATER))
    assert cb.fluorescence_enhancement([1, 1, 1], [R, R, R], 1.0) > 0
    # Strahlen: Leistung 1, Poynting-Fluss 1 (wie test_beams)
    beam = cb.BeamField.gaussian(WATER, 2.0, (0, 0, 0), 3.0, (1, 0, 0), 24, 48)
    assert abs(beam.power - 1) < 1e-12
    assert abs(beam.poynting_flux(0.5, 12.0, 120) - 1) < 2e-3
    E0, H0 = beam.eval((0, 0, 0))
    assert np.isclose(np.sum(np.abs(E0) ** 2), beam.reference_E2)
    ev = cb.near_field_evaluator(m, r.h, b, WATER, om, cb.DipoleField(WATER, om, x0, p))
    assert ev(np.array([[0, 0, 3.0]]))["E"].shape == (1, 3)


# --- Geschichtete Probleme ------------------------------------------------------------------------------------------------------
def test_layered_neutral():
    """Neutrale Schicht (Medium = Aussenmedium): Duennschicht exakt ohne Wirkung, Zweitor und exakte Rechnung bis auf die
    Diskretisierung der zusaetzlichen Flaeche (wie tests/test_thin_layer.cpp und tests/test_twoport.cpp)."""
    om, m, so = 0.5, cb.make_icosphere(4), cb.SolveOptions(tol=1e-9)
    bare = cb.ScatteringProblem(m, GOLD, om, outer=WATER).solve_plane_wave(Z, X, so).sigma_ext
    c = [cb.Coating(0.05, WATER)]
    s_thin = cb.ThinLayerScatteringProblem(m, GOLD, c, om, outer=WATER).solve_plane_wave(Z, X, so).sigma_ext
    assert abs(s_thin / bare - 1) < 1e-10
    body = cb.ThinBody(m, GOLD, c)
    assert len(body.coatings) == 1 and body.inner_fraction == 0.0
    assert np.isclose(cb.ThinLayerScatteringProblem([body], om, outer=WATER).solve_plane_wave(Z, X, so).sigma_ext, s_thin)
    g = cb.LayeredGeometry(WATER)
    regions = g.add_coated_body(m, GOLD, c)
    assert len(regions) == 2 and len(g.surfaces) == 2 and g.triangles == 640
    s_exact = cb.LayeredScatteringProblem(g, om).solve_plane_wave(Z, X, so).sigma_ext
    assert abs(s_exact / bare - 1) < 0.02, s_exact / bare - 1
    # Zweitor: duenne neutrale Schicht d = 1e-3 wie ohne Schicht
    c = [cb.Coating(1e-3, WATER)]
    surf = cb.TwoPortLayerProblem.layer_surfaces(m, c)
    tp = cb.TwoPortLayerProblem(surf, GOLD, c, om, outer=WATER)
    s_tp = tp.solve_plane_wave(Z, X, so).sigma_ext
    assert len(surf) == 2 and tp.layers() >= 1 and abs(s_tp / bare - 1) < 1e-3, s_tp / bare - 1
    tb = cb.TwoPortBody(surf, GOLD, c)
    assert len(tb.layers) == 1 and raises(ValueError, cb.TwoPortBody, surf, GOLD, [])
    assert np.isclose(cb.TwoPortLayerProblem([tb], om, outer=WATER).solve_plane_wave(Z, X, so).sigma_ext, s_tp)


# --- Fehlerbehandlung, GMRES, Threads ---------------------------------------------------------------------------------------------
def test_errors_and_gmres():
    m = cb.make_icosphere(3)
    P = cb.ScatteringProblem(m, GOLD, 0.5)
    assert raises(ValueError, P.solve_rhs, np.zeros(10))
    assert raises(ValueError, cb.exterior_near_field, m, np.zeros(7), WATER, 0.5, Z, X, np.zeros((1, 3)))
    assert raises(ValueError, cb.exterior_near_field, m, np.zeros(8 * len(m)), WATER, 0.5, Z, X, np.zeros((4, 2)))
    assert raises(TypeError, P.solve_plane_wave, "z", X)
    # GMRES mit Python-Operator und Python-Vorkonditionierer
    d = np.linspace(1, 3, 40) + 0.1j
    rhs = np.ones(40, dtype=complex)
    x, info = cb.gmres(lambda v: d * v, rhs, M=lambda v: v / d, tol=1e-12)
    assert info.converged and np.allclose(d * x, rhs)
    assert raises(ValueError, cb.gmres, lambda v: v[:3], rhs)            # falsche Ergebnislaenge
    def boom(v):
        raise ZeroDivisionError("aus Python")
    assert raises(ZeroDivisionError, cb.gmres, boom, rhs)                 # Ausnahmen kommen durch


def test_threads_release_gil():
    """Zwei Loesungen parallel in Python-Threads (ohne GIL) liefern dasselbe wie nacheinander."""
    m = cb.make_icosphere(4)
    P = cb.ScatteringProblem(m, GOLD, 0.5)
    ref = [P.solve_plane_wave(Z, p).sigma_ext for p in (X, (0, 1, 0))]
    out = [None, None]
    ticks = []
    def work(i, p):
        out[i] = P.solve_plane_wave(Z, p).sigma_ext
    def ticker():                                                         # laeuft nur, wenn der GIL frei ist
        t0 = time.time()
        while any(o is None for o in out) and time.time() - t0 < 60:
            ticks.append(1)
            time.sleep(0.001)
    th = [threading.Thread(target=work, args=(0, X)), threading.Thread(target=work, args=(1, (0, 1, 0))), threading.Thread(target=ticker)]
    for t in th:
        t.start()
    for t in th:
        t.join()
    assert np.allclose(out, ref, rtol=1e-12) and len(ticks) > 5


# --- eigene einfallende Felder (v0.43) -----------------------------------------------------------------------------------------
class CountingField(cb.CustomField):
    """Python-Kopie einer Feldklasse mit Zaehler der fields-Aufrufe (ein Aufruf je Projektion bzw. Nahfeldauswertung)."""

    def __init__(self, inner):
        self.inner, self.calls = inner, 0
        self.reference_E2, self.reference_C = inner.reference_E2, inner.reference_C

    def fields(self, x):
        self.calls += 1
        return self.inner.fields(x)


class PyDipole(cb.CustomField):
    """Dipolfeld nach der Formel in dipole.hpp, unabhaengig vom Kern in NumPy."""

    def __init__(self, medium, omega, r0, p):
        self.k, self.eps, self.om = medium.k(omega), medium.eps, omega
        self.r0, self.p = np.asarray(r0, float), np.asarray(p, complex)

    def fields(self, x):
        R = x - self.r0
        r = np.linalg.norm(R, axis=1)[:, None]
        n = R / r
        k = self.k
        G = np.exp(1j * k * r) / (4 * np.pi * r)
        nxp = np.cross(n, self.p)
        H = self.om * k * nxp * G * (1 - 1 / (1j * k * r))
        npn = np.sum(n * self.p, axis=1)[:, None]
        E = (k * k * np.cross(nxp, n) + (3 * n * npn - self.p) * (1 / r ** 2 - 1j * k / r)) * G / self.eps
        return E, H


def test_custom_field_identical_to_core():
    """Python-Welle = PlaneWaveField des Kerns: Projektion, Nahfeld, Kraft und Gradienten bitgleich; ein fields-Aufruf je Stapel."""
    om, m = 0.5, cb.make_icosphere(4)
    core = cb.PlaneWaveField(WATER, om, Z, X)
    py = CountingField(cb.PythonPlaneWave(WATER, om, Z, X))
    assert py.reference_E2 == core.reference_E2 and np.isclose(py.reference_C, core.reference_C)
    b = py.project(m, WATER)
    assert py.calls == 1 and np.array_equal(b, core.project(m, WATER))
    assert np.allclose(b, cb.plane_wave_trace(m, WATER, om, Z, X), atol=1e-14)
    P = cb.ScatteringProblem(m, GOLD, om, outer=WATER)
    r = P.solve_rhs(b, cb.SolveOptions(tol=1e-10))
    pts = np.array([[1.5, 0, 0], [0, 0.3, 2.0], [0.0, 0.0, 0.0]])
    n_core = cb.exterior_near_field(P.mesh, r.h, b, WATER, om, core, pts)
    n_py = cb.exterior_near_field(P.mesh, r.h, b, WATER, om, py, pts)
    assert py.calls == 2 and list(n_py["inside"]) == [False, False, True]
    for key in ("E", "H", "enhancement", "chirality"):
        assert np.array_equal(n_core[key], n_py[key]), key
    ev_core = cb.near_field_evaluator(P.mesh, r.h, b, WATER, om, core)
    ev_py = cb.near_field_evaluator(P.mesh, r.h, b, WATER, om, py)
    del P, r, b                                                       # Auswerter haelt Netz, Spuren und Feld selbst
    py.calls = 0
    F_py = cb.force_on_sphere(ev_py, WATER, (0, 0, 0), 1.5)
    assert py.calls == 1 and np.array_equal(F_py, cb.force_on_sphere(ev_core, WATER, (0, 0, 0), 1.5))
    g_py = cb.fields_with_gradients(ev_py, np.array([[2.0, 0, 0]]), 1e-3)[0]
    g_core = cb.fields_with_gradients(ev_core, np.array([[2.0, 0, 0]]), 1e-3)[0]
    assert py.calls == 2 and np.array_equal(g_py.dE, g_core.dE)
    assert np.array_equal(cb.force_on_offset(ev_py, WATER, m, 0.2), cb.force_on_offset(ev_core, WATER, m, 0.2))


def test_custom_dipole_against_core():
    om, x0, p = 0.5, (0.2, -0.1, 1.4), (0.3, 1j, 1.0)
    f = PyDipole(WATER, om, x0, p)
    x = np.random.default_rng(2).normal(size=(40, 3)) * 2
    E, H = f.eval(x)
    Ec, Hc = cb.dipole_field(WATER, om, x0, p, x)
    assert np.allclose(E, Ec, rtol=1e-12, atol=0) and np.allclose(H, Hc, rtol=1e-12, atol=0)
    m = cb.make_icosphere(4)
    assert np.allclose(f.project(m, WATER), cb.DipoleField(WATER, om, x0, p).project(m, WATER), rtol=1e-12, atol=1e-15)
    # Nahfeld des Dipols variiert auf der Skala r < Wellenlaenge: kleinerer Differenzenschritt (Fehler ~ (h/r)^4)
    assert cb.maxwell_residual(f, WATER, om, x) < 1e-5 and cb.maxwell_residual(f, WATER, om, x, step=2e-3) < 1e-8


def test_superposition():
    """Linearitaet: Loesung und Nahfeld der Summe = Summe; Stehwelle: keine Kraft auf die Kugel im Bauch (Symmetrie)."""
    om, m = 0.5, cb.make_icosphere(4)
    so = cb.SolveOptions(tol=1e-11)
    up, down = cb.PlaneWaveField(WATER, om, Z, X), cb.PlaneWaveField(WATER, om, (0, 0, -1), X)
    dip = PyDipole(WATER, om, (0, 0, 2.5), (1, 0, 0))
    mix = up + 0.5j * dip - down                                        # Kern- und Python-Felder gemischt
    assert isinstance(mix, cb.SuperposedField) and len(mix.components) == 3 and mix.coefficients == [1, 0.5j, -1]
    assert len((mix + up).components) == 4 and isinstance(sum([up, down]), cb.SuperposedField)
    P = cb.ScatteringProblem(m, GOLD, om, outer=WATER)
    parts = [up, dip, down]
    bs = [f.project(m, WATER) for f in parts]
    hs = [P.solve_rhs(b, so).h for b in bs]
    b = mix.project(m, WATER)
    assert np.allclose(b, bs[0] + 0.5j * bs[1] - bs[2], atol=1e-14)
    h = P.solve_rhs(b, so).h
    assert np.allclose(h, hs[0] + 0.5j * hs[1] - hs[2], atol=1e-8 * np.abs(h).max())
    pts = np.array([[1.4, 0.2, 0.1], [0.0, -1.6, 0.4]])
    nf = cb.exterior_near_field(P.mesh, h, b, WATER, om, mix, pts)["E"]
    ref = sum(c * cb.exterior_near_field(P.mesh, hh, bb, WATER, om, f, pts)["E"] for c, f, hh, bb in zip([1, 0.5j, -1], parts, hs, bs))
    assert np.allclose(nf, ref, atol=1e-8)
    # Stehwelle E = 2 x cos(kz): spiegelsymmetrisch zu z = 0, also F_z = 0 (Rauschboden der Quadratur: einige 1e-6 der
    # Einzelwelle, so gross wie deren Querkraft); um +-s verschoben: F_z(s) = -F_z(-s), deutlich von null verschieden
    def standing(s):
        sw = cb.PythonPlaneWave(WATER, om, Z, X, origin=(0, 0, s)) + cb.PythonPlaneWave(WATER, om, (0, 0, -1), X, origin=(0, 0, s))
        b = sw.project(m, WATER)
        r = P.solve_rhs(b, so)
        return cb.force_on_sphere(cb.near_field_evaluator(P.mesh, r.h, b, WATER, om, sw), WATER, (0, 0, 0), 1.5)
    F1 = cb.force_on_sphere(P.mesh, hs[0], WATER, om, Z, X, (0, 0, 0), 1.5)
    assert F1[2] > 0 and abs(F1[0]) < 2e-5 * F1[2]
    assert abs(standing(0.0)[2]) < 2e-5 * F1[2]
    Fp, Fm = standing(0.3), standing(-0.3)
    assert abs(Fp[2]) > 0.05 * F1[2] and abs(Fp[2] + Fm[2]) < 1e-4 * abs(Fp[2]), (Fp, Fm)


def test_maxwell_residual_and_beams():
    om = 0.5
    pts = np.random.default_rng(3).normal(size=(20, 3))
    assert cb.maxwell_residual(cb.PythonPlaneWave(WATER, om, (1, 1, 1), (1, -1, 0)), WATER, om, pts) < 1e-8
    chiral = cb.Medium(eps=1.7689, chi=0.02)
    beam = cb.BeamField.gaussian(chiral, 2.0, (0, 0, 0), 3.0, cb.circular_polarization(Z, +1), 24, 48)
    assert cb.maxwell_residual(beam, chiral, 2.0, pts) < 1e-8          # Pasteur-Relationen D, B wie im Kern
    pw = cb.PythonPlaneWave(WATER, om, Z, X)
    wrong_H = cb.as_field_function(lambda x: (pw.fields(x)[0], 2 * pw.fields(x)[1]))
    wrong_k = cb.PythonPlaneWave(cb.Medium(eps=2.25), om, Z, X)
    assert cb.maxwell_residual(wrong_H, WATER, om, pts) > 0.1 and cb.maxwell_residual(wrong_k, WATER, om, pts) > 0.05
    # Bessel-Strahl: exakte Loesung, Intensitaet J0^2 auf der Achse unabhaengig von z; Wirbel (Ladung 1) dunkel auf der Achse
    bes = cb.AngularSpectrumField.bessel(WATER, om, 0.6, nphi=48)
    assert cb.maxwell_residual(bes, WATER, om, pts) < 1e-8
    ax = np.array([[0, 0, z] for z in (-3.0, 0.0, 2.0)])
    I_ax = np.sum(np.abs(bes.eval(ax)[0]) ** 2, axis=1)
    assert np.allclose(I_ax, I_ax[0]) and np.isclose(I_ax[0], bes.reference_E2)
    # Wirbel, Ladung 1: auf der Achse ueberlebt nur die Harmonische 0 der Amplituden. Zirkular mit sigma = +1: dunkel;
    # sigma = -1: Spin-Bahn-Kopplung, nur E_z auf der Achse (linear waere ebenfalls nur E_z, nicht dunkel)
    def vortex(sg):
        return cb.AngularSpectrumField.bessel(WATER, om, 0.6, pol=cb.circular_polarization(Z, sg) / np.sqrt(2), charge=1, nphi=48)
    off = np.sum(np.abs(vortex(+1).eval((1.5, 0, 1.0))[0]) ** 2)
    assert np.sum(np.abs(vortex(+1).eval((0, 0, 1.0))[0]) ** 2) < 1e-20 * off
    E_ax = vortex(-1).eval((0, 0, 1.0))[0]
    assert np.abs(E_ax[:2]).max() < 1e-12 and abs(E_ax[2]) > 0.1
    assert raises(ValueError, cb.AngularSpectrumField, chiral, om, [Z], [X])


def test_custom_field_errors():
    om, m = 0.5, cb.make_icosphere(3)
    P = cb.ScatteringProblem(m, GOLD, om, outer=WATER)
    pw = cb.PythonPlaneWave(WATER, om, Z, X)
    b = pw.project(m, WATER)
    r = P.solve_rhs(b)
    pts = np.array([[2.0, 0, 0]])
    bad_shape = cb.as_field_function(lambda x: (np.zeros((len(x), 2)), np.zeros((len(x), 3))))
    not_pair = cb.as_field_function(lambda x: np.zeros((len(x), 3)))
    assert raises(ValueError, cb.exterior_near_field, m, r.h, b, WATER, om, bad_shape, pts)
    assert raises(TypeError, cb.exterior_near_field, m, r.h, b, WATER, om, not_pair, pts)
    assert raises(ValueError, bad_shape.project, m, WATER)
    assert raises(TypeError, cb.exterior_near_field, m, r.h, b, WATER, om, object(), pts)
    assert raises(NotImplementedError, cb.CustomField().project, m, WATER)
    as_list = cb.as_field_function(lambda x: list(pw.fields(x)))     # Liste statt Tupel ist erlaubt
    assert np.array_equal(cb.exterior_near_field(m, r.h, b, WATER, om, as_list, pts)["E"],
                          cb.exterior_near_field(m, r.h, b, WATER, om, pw, pts)["E"])

    class Boom(cb.CustomField):
        def fields(self, x):
            raise ZeroDivisionError("aus Python")
    ev = cb.near_field_evaluator(m, r.h, b, WATER, om, Boom())
    assert raises(ZeroDivisionError, cb.force_on_sphere, ev, WATER, (0, 0, 0), 1.5)   # Ausnahme aus dem Lauf ohne GIL
    assert np.all(np.isfinite(cb.force_on_sphere(cb.near_field_evaluator(m, r.h, b, WATER, om, pw), WATER, (0, 0, 0), 1.5)))


# --- lineare Dichten (v0.45) -----------------------------------------------------------------------------------------------
def test_linear_densities():
    """Unstetig lineare Dichten: Summenidentitaet der Eintraege, Projektion, Streuung wie test_linear (C++) und wie die
    Vorhersage aus Stufe 1 (eben/linear, Glas, 320 Elemente: -6,163 %), Fehlerbehandlung, Lebensdauer."""
    m = cb.make_icosphere(4)
    E, K = cb.LinearKernelEntries(m, 1.3 + 0.05j), cb.KernelEntries(m, 1.3 + 0.05j)
    for i, j in [(0, 0), (0, 1), (0, 50), (3, 200)]:
        assert np.abs(E.lambda_block(i, j).sum(axis=(0, 1)) - K.exact(i, j)).max() < 1e-14 * np.abs(K.exact(i, j)).max()
        assert np.allclose(E.block(i, j).sum(axis=(0, 1)) / 3, K.exact(i, j) / np.sqrt(m.areas[i] * m.areas[j]), rtol=1e-12)
    S = E.S(5)
    lam_gram = m.areas[5] / 12 * (np.eye(3) + 1)
    assert np.allclose(S @ lam_gram @ S.T, np.eye(3))                  # psi orthonormal
    b = cb.plane_wave_trace_linear(m, cb.Medium(), 1.0, Z, X)
    assert np.abs(cb.linear_to_constant(m, b) - cb.plane_wave_trace(m, cb.Medium(), 1.0, Z, X)).max() < 1e-14
    P = cb.LinearScatteringProblem(m, cb.Medium(eps=2.25), 1.0)
    r = P.solve_plane_wave(Z, X, cb.SolveOptions(tol=1e-8))
    err = r.sigma_ext / np.pi / 0.2150978 - 1
    assert P.unknowns == 24 * len(m) == len(r.h) and abs(err - (-0.06163)) < 1e-4, err
    hs = r.h - b
    assert np.isclose(cb.extinction_cross_section_linear(m, hs, 1.0, 1.0, Z, X), r.sigma_ext)
    assert np.isclose(r.sigma_ext, 4 * np.pi * r.forward.real)
    r2 = P.solve_rhs(b, cb.SolveOptions(tol=1e-8))
    assert np.allclose(r2.h, r.h)
    v = cb.linear_trace_value(m, r.h, 7, m.centroids[7])                # Dichte im Schwerpunkt = Mittelwert
    assert np.allclose(v.c, cb.linear_to_constant(m, r.h)[56:64] / np.sqrt(m.areas[7]))
    assert raises(ValueError, P.solve_rhs, np.zeros(8 * len(m)))
    assert raises(ValueError, cb.plane_wave_trace_linear, m, cb.Medium(eps=2.0, chi=0.1), 1.0, Z, cb.circular_polarization(Z, 1))
    H = cb.linear_hmatrix(cb.LinearKernelEntries(m, 1.0))                # Eintraege und Netz bleiben am Leben
    C = cb.LinearCauchyOperator(m, H)
    del H
    # Plemelj E b = b (ebene Welle = innere Loesung): k = 1, 320 Elemente: linear 1,1e-3, konstant 7,3e-3
    Cc = cb.CauchyOperator(m, cb.KernelHMatrix(cb.KernelEntries(m, 1.0)))
    bc = cb.plane_wave_trace(m, cb.Medium(), 1.0, Z, X)
    e_lin, e_con = (np.linalg.norm(C.apply(b) - b) / np.linalg.norm(b), np.linalg.norm(Cc.apply(bc) - bc) / np.linalg.norm(bc))
    assert C.apply(b).shape == b.shape and e_lin < 2e-3 and e_lin < 0.25 * e_con, (e_lin, e_con)


def main():
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_") and callable(f)]
    failed = 0
    for name, f in tests:
        t0 = time.time()
        try:
            f()
            print(f"  ok      {name}  ({time.time() - t0:.1f} s)")
        except Exception:
            failed += 1
            print(f"  FEHLER  {name}")
            traceback.print_exc()
    print(f"cliffordbem {cb.__version__}: {len(tests) - failed}/{len(tests)} Tests bestanden")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
