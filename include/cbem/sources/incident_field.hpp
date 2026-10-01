#pragma once
// Allgemeine einfallende Felder (v0.37): ebene Welle, Dipolquelle, Strahl als Winkelspektrum ebener Wellen (Gaussstrahl,
// fokussierter Strahl nach Richards-Wolf). Jede Ueberlagerung ebener Wellen ist eine exakte Loesung der Maxwell-Gleichungen im
// homogenen Medium -- ohne paraxiale Naeherung, fuer die BEM direkt verwendbar.
//
// Strahl: E(r) = int a(k) e^{ik k.(r - r_f)} dOmega (k Einheitsvektor, Brennpunkt r_f, Achse +z), H = sqrt(eps/mu) k x E je Teilwelle.
//   Richards-Wolf: a = f(theta) sqrt(cos theta) [(E_in.phi) phi + (E_in.rho) theta_hat] fuer theta <= theta_max = asin(NA/n),
//                  f(theta) = exp(-sin^2 theta / (f0^2 sin^2 theta_max)) (gaussfoermige Pupillenausleuchtung, Fuellfaktor f0)
//   Gaussstrahl:   a = exp(-(k w0 sin theta / 2)^2) (p_x cos theta, p_y cos theta, -sin theta (p_x cos phi + p_y sin phi))
//                  (exakt aus dem Querfeld exp(-rho^2/w0^2) p im Brennpunkt, jede Teilwelle transversal)
// Leistung (Parseval): P = 1/2 sqrt(eps/mu) (2 pi/k)^2 int |a|^2 dOmega; der Strahl wird auf P = 1 normiert.
// Chirales Medium (v0.38): jede Teilwelle a(k) wird in ihre Helizitaetsanteile c_s e_s(k) zerlegt (e_s = circular_polarization/
// sqrt 2, orthonormal), jeder Anteil laeuft mit seiner Wellenzahl k_s; H = sqrt(eps/mu) k x E je Anteil (die Wellenimpedanz ist im
// Pasteur-Medium unveraendert); P = sum_s 1/2 sqrt(eps/mu) (2 pi/k_s)^2 int |c_s|^2 dOmega. Kraft-Effizienz der
// optischen Pinzette Q = F c / (n P) = F / sqrt(eps mu) (c = 1).
#include <memory>
#include <vector>
#include "cbem/sources/chiral_incidence.hpp"

namespace cbem {

class IncidentField {
public:
    virtual ~IncidentField() = default;
    virtual void eval(const Vec3& x, CVec3& E, CVec3& H) const = 0;
    virtual real reference_E2() const = 0;      // Bezug fuer die Feldverstaerkung |E|^2 / |E0|^2
    virtual real reference_C() const = 0;       // Bezug fuer die optische Chiralitaet (|Im(E0*.H0)| einer zirkularen Welle)
    // Spur F = sqrt(eps) E + I sqrt(mu) H auf den stueckweise konstanten Dichten (7-Punkt-Regel, sub-fach unterteilt)
    std::vector<cplx> project(const TriangleMesh& mesh, const Medium& m, int sub = 2) const;
};

class PlaneWaveField : public IncidentField {
public:
    PlaneWaveField(const Medium& m, real omega, const Vec3& d, const CVec3& p);
    void eval(const Vec3& x, CVec3& E, CVec3& H) const override;
    real reference_E2() const override { return p2_; }
    real reference_C() const override { return C0_; }
private:
    Medium m_; Vec3 d_; CVec3 p_, dxp_; cplx k_; real p2_, C0_;
};

class ZeroField : public IncidentField {                  // kein einfallendes Feld: das Nahfeld liefert dann das Streufeld allein
public:
    void eval(const Vec3&, CVec3& E, CVec3& H) const override { E = CVec3{}; H = CVec3{}; }
    real reference_E2() const override { return 1.0; }
    real reference_C() const override { return 1.0; }
};

class DipoleField : public IncidentField {                // elektrischer und magnetischer Dipol (dipole.hpp)
public:
    DipoleField(const Medium& m, real omega, const Vec3& r0, const CVec3& p, const CVec3& md = CVec3{});
    void eval(const Vec3& x, CVec3& E, CVec3& H) const override;
    real reference_E2() const override { return 1.0; }
    real reference_C() const override { return 1.0; }
    const Vec3& position() const { return r0_; }
    const CVec3& p() const { return p_; }
    const CVec3& md() const { return md_; }
private:
    Medium m_; real omega_; Vec3 r0_; CVec3 p_, md_;
};

class BeamField : public IncidentField {                  // Winkelspektrum, auf Leistung 1 normiert
public:
    // fokussierter Strahl (Richards-Wolf): NA, Fuellfaktor f0, Polarisation in der Pupille (p_x, p_y)
    static BeamField focused(const Medium& m, real omega, const Vec3& focus, real NA, real f0, const CVec3& pupil_pol, int ntheta = 40, int nphi = 80);
    // Gaussstrahl mit Taille w0 im Brennpunkt
    static BeamField gaussian(const Medium& m, real omega, const Vec3& focus, real w0, const CVec3& pol, int ntheta = 40, int nphi = 80);
    void eval(const Vec3& x, CVec3& E, CVec3& H) const override;
    real reference_E2() const override { return I0_; }    // |E|^2 im Brennpunkt
    real reference_C() const override { return I0_ * std::real(std::sqrt(m_.eps / m_.mu)); }
    real power() const { return P_; }                     // nach der Normierung 1
    std::size_t components() const { return dir_.size(); }
    // unabhaengige Pruefung: Poynting-Fluss durch die Ebene z = z_f + dz, Quadrat +-half_width, ng x ng Mittelpunktsregel
    real poynting_flux(real dz, real half_width, int ng) const;
private:
    BeamField(const Medium& m, real omega, const Vec3& focus) : m_(m), omega_(omega), focus_(focus) {}
    void finish();                                        // Leistung, Normierung, Bezugsintensitaet
    Medium m_; real omega_; Vec3 focus_; cplx k_ = 0;
    std::vector<Vec3> dir_; std::vector<CVec3> amp_;      // a(k) dOmega je Teilwelle
    std::vector<real> dw_;                                // dOmega je Teilwelle (fuer die Leistung nach Parseval)
    std::vector<cplx> kc_;                                // Wellenzahl je Teilwelle (chiral: k_sigma der Helizitaet)
    void add(const Vec3& dir, const CVec3& a, real w);    // Teilwelle anfuegen; im chiralen Medium nach Helizitaeten zerlegt
    real P_ = 0, I0_ = 0;
};

}  // namespace cbem
