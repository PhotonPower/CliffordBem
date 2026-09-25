#include "cbem/core/materials.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cbem {

std::string ConstantMaterial::name() const {
    std::ostringstream s; s << "konstant(" << e_.real() << (e_.imag() >= 0 ? "+" : "") << e_.imag() << "i)"; return s.str();
}

TabulatedMaterial::TabulatedMaterial(std::string name, std::vector<real> l, std::vector<real> n, std::vector<real> k)
    : name_(std::move(name)), l_(std::move(l)), n_(std::move(n)), k_(std::move(k)) {
    if (l_.size() < 2 || l_.size() != n_.size() || l_.size() != k_.size()) throw std::invalid_argument("TabulatedMaterial: Daten");
}

std::shared_ptr<TabulatedMaterial> TabulatedMaterial::from_yaml(const std::string& path) {
    std::ifstream in(path); if (!in) throw std::runtime_error("Materialdatei nicht lesbar: " + path);
    std::string line; bool data = false; std::vector<real> l, n, k;
    while (std::getline(in, line)) {
        if (line.find("tabulated nk") != std::string::npos) { data = false; continue; }
        if (line.find("data: |") != std::string::npos) { data = true; continue; }
        if (!data) continue;
        std::istringstream s(line); real a, b, c;
        if (s >> a >> b >> c) { l.push_back(1000.0 * a); n.push_back(b); k.push_back(c); }
        else if (!line.empty() && line.find_first_not_of(" \t") != std::string::npos) data = false;
    }
    std::string nm = path.substr(path.find_last_of('/') + 1);
    return std::make_shared<TabulatedMaterial>(nm, l, n, k);
}

cplx TabulatedMaterial::eps(real lam) const {
    if (lam < l_.front() || lam > l_.back()) throw std::out_of_range("Wellenlaenge ausserhalb der Tabelle von " + name_);
    std::size_t i = std::upper_bound(l_.begin(), l_.end(), lam) - l_.begin();
    if (i == 0) i = 1;
    if (i >= l_.size()) i = l_.size() - 1;
    real t = (lam - l_[i - 1]) / (l_[i] - l_[i - 1]);
    cplx nk(n_[i - 1] + t * (n_[i] - n_[i - 1]), k_[i - 1] + t * (k_[i] - k_[i - 1]));
    return nk * nk;
}

std::shared_ptr<Material> make_material(const std::string& spec, const std::string& dir) {
    if (spec == "Au" || spec == "Ag") return TabulatedMaterial::from_yaml(dir + "/" + spec + "_Johnson.yml");
    auto c = spec.find(',');
    return std::make_shared<ConstantMaterial>(cplx(std::stod(spec.substr(0, c)), c == std::string::npos ? 0.0 : std::stod(spec.substr(c + 1))));
}

}  // namespace cbem
