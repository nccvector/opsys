#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <variant>

namespace osys {

struct ConstantIndex {
    double n{};
};

struct Sellmeier3 {
    // C coefficients are in micrometer^2. Input wavelengths are still nm.
    std::array<double, 3> b{};
    std::array<double, 3> c_um2{};
};

struct Medium {
    constexpr Medium() = default;
    explicit constexpr Medium(ConstantIndex constant) : model(constant) {}
    explicit constexpr Medium(Sellmeier3 sellmeier) : model(sellmeier) {}

    [[nodiscard]] double refractive_index(double wavelength_nm) const;

    std::variant<ConstantIndex, Sellmeier3> model{ConstantIndex{1.0}};
};

[[nodiscard]] inline double refractive_index(const Medium& medium, double wavelength_nm) {
    return std::visit([wavelength_nm](const auto& model) {
        using Model = std::decay_t<decltype(model)>;

        if constexpr (std::is_same_v<Model, ConstantIndex>) {
            return model.n;
        } else {
            const double lambda_um = wavelength_nm * 0.001;
            const double lambda2 = lambda_um * lambda_um;
            double n2 = 1.0;

            for (std::size_t i = 0; i < model.b.size(); ++i) {
                n2 += model.b[i] * lambda2 / (lambda2 - model.c_um2[i]);
            }

            return std::sqrt(n2);
        }
    }, medium.model);
}

inline double Medium::refractive_index(double wavelength_nm) const {
    return osys::refractive_index(*this, wavelength_nm);
}

inline constexpr Medium air_medium{ConstantIndex{1.0}};

inline constexpr Medium n_bk7_medium{Sellmeier3{
        .b = {1.03961212, 0.231792344, 1.01046945},
        .c_um2 = {0.00600069867, 0.0200179144, 103.560653},
}};

inline constexpr Medium dense_medium{ConstantIndex{1.5}};

} // namespace osys
