#pragma once
// Gemeinsame Hilfen der Python-Anbindung (v0.42): Typumwandler fuer Vec3/CVec3, NumPy-Konvertierungen, Parameterklassen.
//
// Konventionen:
//  - Vec3 und CVec3 nehmen jede Folge von drei Zahlen an (Tupel, Liste, numpy.ndarray) und werden als numpy.ndarray der Laenge 3
//    zurueckgegeben (float64 bzw. complex128).
//  - Spuren (8 Koeffizienten je Dreieck) sind flache complex128-Arrays der Laenge 8 N wie im C++-Kern; reshape(-1, 8) gibt die
//    Multivektor-Koeffizienten je Dreieck (Blade-Reihenfolge 1, e1, e2, e12, e3, e13, e23, e123).
//  - Punktwolken sind float64-Arrays der Form (M, 3); ein einzelner Punkt darf als Form (3,) uebergeben werden.
//  - Mat8 (Linksmultiplikation, Systembloecke) wird als (8, 8)-Array zurueckgegeben, Zeile r, Spalte c.
// Diese Datei muss vor jeder Verwendung der Typen in einer Uebersetzungseinheit eingebunden werden (Spezialisierungen der
// type_caster haben Vorrang vor dem allgemeinen std::array-Umwandler aus pybind11/stl.h).
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstring>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "cbem/clifford/multivector.hpp"
#include "cbem/core/types.hpp"
#include "cbem/linalg/dense.hpp"

namespace pybind11 {
namespace detail {

template <> struct type_caster<cbem::Vec3> {
    PYBIND11_TYPE_CASTER(cbem::Vec3, const_name("numpy.ndarray[float64[3]]"));
    bool load(handle src, bool) {
        if (!src || isinstance<str>(src) || isinstance<bytes>(src) || !PySequence_Check(src.ptr())) return false;
        auto seq = reinterpret_borrow<sequence>(src);
        if (seq.size() != 3) return false;
        try {
            for (int i = 0; i < 3; ++i) value[i] = seq[static_cast<size_t>(i)].cast<double>();
        } catch (const cast_error&) { return false; }
        return true;
    }
    static handle cast(const cbem::Vec3& v, return_value_policy, handle) {
        array_t<double> a(3);
        double* d = a.mutable_data();
        d[0] = v.x; d[1] = v.y; d[2] = v.z;
        return a.release();
    }
};

template <> struct type_caster<cbem::CVec3> {
    PYBIND11_TYPE_CASTER(cbem::CVec3, const_name("numpy.ndarray[complex128[3]]"));
    bool load(handle src, bool) {
        if (!src || isinstance<str>(src) || isinstance<bytes>(src) || !PySequence_Check(src.ptr())) return false;
        auto seq = reinterpret_borrow<sequence>(src);
        if (seq.size() != 3) return false;
        try {
            for (int i = 0; i < 3; ++i) value[static_cast<size_t>(i)] = seq[static_cast<size_t>(i)].cast<std::complex<double>>();
        } catch (const cast_error&) { return false; }
        return true;
    }
    static handle cast(const cbem::CVec3& v, return_value_policy, handle) {
        array_t<std::complex<double>> a(3);
        std::copy(v.begin(), v.end(), a.mutable_data());
        return a.release();
    }
};

}  // namespace detail
}  // namespace pybind11

namespace cbem::py_bind {

namespace py = pybind11;
using namespace pybind11::literals;

using CArr = py::array_t<cplx, py::array::c_style | py::array::forcecast>;
using RArr = py::array_t<double, py::array::c_style | py::array::forcecast>;

// --- C++ -> NumPy -----------------------------------------------------------------------------------------------------------
inline py::array_t<cplx> to_numpy(const std::vector<cplx>& v) {
    py::array_t<cplx> a(static_cast<py::ssize_t>(v.size()));
    if (!v.empty()) std::memcpy(a.mutable_data(), v.data(), v.size() * sizeof(cplx));
    return a;
}
inline py::array_t<double> to_numpy(const std::vector<real>& v) {
    py::array_t<double> a(static_cast<py::ssize_t>(v.size()));
    if (!v.empty()) std::memcpy(a.mutable_data(), v.data(), v.size() * sizeof(real));
    return a;
}
inline py::array_t<double> points_to_numpy(const std::vector<Vec3>& v) {
    py::array_t<double> a({static_cast<py::ssize_t>(v.size()), py::ssize_t(3)});
    double* d = a.mutable_data();
    for (std::size_t i = 0; i < v.size(); ++i) { d[3 * i] = v[i].x; d[3 * i + 1] = v[i].y; d[3 * i + 2] = v[i].z; }
    return a;
}
inline py::array_t<cplx> cvecs_to_numpy(const std::vector<CVec3>& v) {
    py::array_t<cplx> a({static_cast<py::ssize_t>(v.size()), py::ssize_t(3)});
    cplx* d = a.mutable_data();
    for (std::size_t i = 0; i < v.size(); ++i) for (int c = 0; c < 3; ++c) d[3 * i + c] = v[i][c];
    return a;
}
inline py::array_t<cplx> mat8_to_numpy(const Mat8& M) {
    py::array_t<cplx> a({py::ssize_t(8), py::ssize_t(8)});
    std::copy(M.begin(), M.end(), a.mutable_data());
    return a;
}
inline py::array_t<cplx> matrix_to_numpy(const Matrix& M) {          // Matrix ist spaltenweise gespeichert
    py::array_t<cplx> a({static_cast<py::ssize_t>(M.rows), static_cast<py::ssize_t>(M.cols)});
    cplx* d = a.mutable_data();
    for (std::size_t r = 0; r < M.rows; ++r) for (std::size_t c = 0; c < M.cols; ++c) d[r * M.cols + c] = M(r, c);
    return a;
}
inline py::array_t<cplx> multivectors_to_numpy(const std::vector<Multivector>& v) {
    py::array_t<cplx> a({static_cast<py::ssize_t>(v.size()), py::ssize_t(8)});
    cplx* d = a.mutable_data();
    for (std::size_t i = 0; i < v.size(); ++i) std::copy(v[i].c.begin(), v[i].c.end(), d + 8 * i);
    return a;
}

// --- NumPy -> C++ -----------------------------------------------------------------------------------------------------------
inline std::vector<cplx> to_cvec(const CArr& a) { return std::vector<cplx>(a.data(), a.data() + a.size()); }
inline std::vector<real> to_rvec(const RArr& a) { return std::vector<real>(a.data(), a.data() + a.size()); }

// Trace der Laenge 8 N pruefen (Spur auf einem Netz mit N Dreiecken)
inline std::vector<cplx> to_trace(const CArr& a, std::size_t N, const char* what) {
    if (static_cast<std::size_t>(a.size()) != 8 * N)
        throw py::value_error(std::string(what) + ": erwartet 8 N = " + std::to_string(8 * N) + " Koeffizienten, erhalten " +
                              std::to_string(a.size()));
    return to_cvec(a);
}

// Punkte (M, 3) oder ein Punkt (3,)
inline std::vector<Vec3> to_points(const RArr& a) {
    if (a.ndim() == 1 && a.shape(0) == 3) return {Vec3(a.data()[0], a.data()[1], a.data()[2])};
    if (a.ndim() != 2 || a.shape(1) != 3) throw py::value_error("Punkte: erwartet ein Array der Form (M, 3) oder (3,)");
    std::vector<Vec3> p(static_cast<std::size_t>(a.shape(0)));
    const double* d = a.data();
    for (std::size_t i = 0; i < p.size(); ++i) p[i] = Vec3(d[3 * i], d[3 * i + 1], d[3 * i + 2]);
    return p;
}
inline bool single_point(const RArr& a) { return a.ndim() == 1; }

// --- GIL ----------------------------------------------------------------------------------------------------------------------
// Rechnung ohne GIL (alle Argumente muessen bereits nach C++ umgewandelt sein)
template <class F> auto nogil(F&& f) -> decltype(f()) {
    py::gil_scoped_release release;
    return f();
}

// --- Parameterklassen ------------------------------------------------------------------------------------------------------
// Bindet eine Struktur mit oeffentlichen Feldern: Konstruktor mit Schluesselwortargumenten (unbekannte Namen -> AttributeError),
// __repr__ mit allen Feldern, Kopie, Gleichheit der Darstellung.
template <class T> class ParamClass {
public:
    ParamClass(py::module_& m, const char* name, const char* doc) : cls_(m, name, doc), names_(std::make_shared<std::vector<std::string>>()) {
        cls_.def(py::init([](py::kwargs kw) {
                     py::object o = py::cast(T{});
                     for (auto item : kw) py::setattr(o, item.first, item.second);
                     return o.cast<T>();
                 }),
                 "Alle Felder sind optional und werden als Schluesselwortargumente gesetzt.");
        auto names = names_;
        const std::string cname = name;
        cls_.def("__repr__", [names, cname](py::object self) {
            std::string s = cname + "(";
            for (std::size_t i = 0; i < names->size(); ++i) {
                if (i) s += ", ";
                s += (*names)[i] + "=" + py::repr(self.attr((*names)[i].c_str())).cast<std::string>();
            }
            return s + ")";
        });
        cls_.def("__copy__", [](const T& t) { return T(t); });
        cls_.def("__deepcopy__", [](const T& t, py::dict) { return T(t); }, "memo"_a);
        cls_.def("copy", [](const T& t) { return T(t); }, "Kopie");
    }
    template <class M> ParamClass& field(const char* n, M T::*p, const char* doc) {
        cls_.def_readwrite(n, p, doc);
        names_->push_back(n);
        return *this;
    }
    py::class_<T>& cls() { return cls_; }
private:
    py::class_<T> cls_;
    std::shared_ptr<std::vector<std::string>> names_;
};

// Teilmodule
void init_core(py::module_& m);
void init_geometry(py::module_& m);
void init_operators(py::module_& m);
void init_problems(py::module_& m);
void init_sources(py::module_& m);
void init_linear(py::module_& m);

}  // namespace cbem::py_bind
