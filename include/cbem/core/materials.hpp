#pragma once
// Materialmodelle: konstante Permittivitaet oder tabellierte optische Konstanten (n, k) ueber der
// Vakuumwellenlaenge, z. B. aus der refractiveindex.info-Datenbank (YAML, "tabulated nk", Wellenlaenge in um).
// eps(lambda) = (n + i k)^2  (Zeitkonvention e^{-i omega t}, Absorption Im eps > 0).
#include <memory>
#include <string>
#include <vector>
#include "cbem/core/types.hpp"

namespace cbem {

class Material {
public:
    virtual ~Material() = default;
    virtual cplx eps(real lambda_nm) const = 0;
    virtual std::string name() const = 0;
};

class ConstantMaterial : public Material {
public:
    explicit ConstantMaterial(cplx e) : e_(e) {}
    cplx eps(real) const override { return e_; }
    std::string name() const override;
private:
    cplx e_;
};

class TabulatedMaterial : public Material {
public:
    // lambda in nm, n, k; lineare Interpolation in n und k
    TabulatedMaterial(std::string name, std::vector<real> lambda_nm, std::vector<real> n, std::vector<real> k);
    static std::shared_ptr<TabulatedMaterial> from_yaml(const std::string& path);   // refractiveindex.info
    cplx eps(real lambda_nm) const override;
    std::string name() const override { return name_; }
    real lambda_min() const { return l_.front(); }
    real lambda_max() const { return l_.back(); }
private:
    std::string name_;
    std::vector<real> l_, n_, k_;
};

// "Au" / "Ag" -> data/materials/<X>_Johnson.yml relativ zu data_dir; sonst komplexe Zahl "re,im"
std::shared_ptr<Material> make_material(const std::string& spec, const std::string& data_dir = "data/materials");

}  // namespace cbem
