#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <variant>

namespace lenses {

struct ConstantIndex {
    double n{};
};

struct Sellmeier3 {
    // C coefficients are in micrometer^2. Input wavelengths are still nm.
    std::array<double, 3> b{};
    std::array<double, 3> c_um2{};
};

struct Material {
    constexpr Material() = default;
    explicit constexpr Material(ConstantIndex constant) : model(constant) {}
    explicit constexpr Material(Sellmeier3 sellmeier) : model(sellmeier) {}

    [[nodiscard]] double refractive_index(double wavelength_nm) const;

    std::variant<ConstantIndex, Sellmeier3> model{ConstantIndex{1.0}};
};

[[nodiscard]] inline double refractive_index(const Material& material, double wavelength_nm) {
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
    }, material.model);
}

inline double Material::refractive_index(double wavelength_nm) const {
    return lenses::refractive_index(*this, wavelength_nm);
}

struct Medium {
    Material material{};
};

[[nodiscard]] inline double refractive_index(const Medium& medium, double wavelength_nm) {
    return refractive_index(medium.material, wavelength_nm);
}

inline constexpr Material air_material{ConstantIndex{1.0}};

inline constexpr Material n_bk7_material{Sellmeier3{
        .b = {1.03961212, 0.231792344, 1.01046945},
        .c_um2 = {0.00600069867, 0.0200179144, 103.560653},
}};

inline constexpr Material dense_material{ConstantIndex{1.5}};

struct FixedMediumCatalog {
    Medium air{};
    Medium n_bk7{};
    Medium dense{};
};

inline constexpr FixedMediumCatalog fixed_medium_catalog{
    .air = Medium{.material = air_material},
    .n_bk7 = Medium{.material = n_bk7_material},
    .dense = Medium{.material = dense_material},
};

} // namespace lenses
