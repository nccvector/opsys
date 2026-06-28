#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <string_view>
#include <type_traits>
#include <variant>

namespace opsys {

struct ConstantIndex {
    double n{};
};

struct Sellmeier3 {
    // C coefficients are in micrometer^2. Input wavelengths are still nm.
    std::array<double, 3> b{};
    std::array<double, 3> c_um2{};
};

struct AbbeVd {
    double nd{};
    double vd{};
};

struct Medium {
    constexpr Medium() = default;
    explicit constexpr Medium(ConstantIndex constant) : model(constant) {}
    explicit constexpr Medium(Sellmeier3 sellmeier) : model(sellmeier) {}
    explicit constexpr Medium(AbbeVd abbe) : model(abbe) {}

    [[nodiscard]] double refractive_index(double wavelength_nm) const;

    std::variant<ConstantIndex, Sellmeier3, AbbeVd> model{ConstantIndex{1.0}};
};

[[nodiscard]] inline double refractive_index(const Medium& medium, double wavelength_nm) {
    return std::visit([wavelength_nm](const auto& model) {
        using Model = std::decay_t<decltype(model)>;

        if constexpr (std::is_same_v<Model, ConstantIndex>) {
            return model.n;
        } else if constexpr (std::is_same_v<Model, Sellmeier3>) {
            const double lambda_um = wavelength_nm * 0.001;
            const double lambda2 = lambda_um * lambda_um;
            double n2 = 1.0;

            for (std::size_t i = 0; i < model.b.size(); ++i) {
                n2 += model.b[i] * lambda2 / (lambda2 - model.c_um2[i]);
            }

            return std::sqrt(n2);
        } else {
            if (model.vd <= 0.0) {
                return model.nd;
            }

            constexpr double lambda_c_um = 0.6562725;
            constexpr double lambda_d_um = 0.5875618;
            constexpr double lambda_f_um = 0.4861327;
            const double lambda_um = wavelength_nm * 0.001;
            const double delta_fc = (model.nd - 1.0) / model.vd;
            const double cauchy_b = delta_fc
                / (1.0 / (lambda_f_um * lambda_f_um) - 1.0 / (lambda_c_um * lambda_c_um));
            const double cauchy_a = model.nd - cauchy_b / (lambda_d_um * lambda_d_um);

            return cauchy_a + cauchy_b / (lambda_um * lambda_um);
        }
    }, medium.model);
}

inline double Medium::refractive_index(double wavelength_nm) const {
    return opsys::refractive_index(*this, wavelength_nm);
}

enum class MediumId {
    air,
    n_bk7,
    n_sf6,
    n_lak9,
    n_sf11,
    dense,
};

inline constexpr Medium air_medium{ConstantIndex{1.0}};

inline constexpr Medium n_bk7_medium{Sellmeier3{
        .b = {1.03961212, 0.231792344, 1.01046945},
        .c_um2 = {0.00600069867, 0.0200179144, 103.560653},
}};

inline constexpr Medium n_sf6_medium{Sellmeier3{
        .b = {1.72448482, 0.390104889, 1.04572858},
        .c_um2 = {0.01349871764, 0.05694099484, 118.5571585},
}};

inline constexpr Medium n_lak9_medium{Sellmeier3{
        .b = {1.46231905, 0.344399589, 1.15508372},
        .c_um2 = {0.00724270156, 0.0243353131, 85.4686868},
}};

inline constexpr Medium n_sf11_medium{Sellmeier3{
        .b = {1.73759695, 0.313747346, 1.89878101},
        .c_um2 = {0.01318871330, 0.06238975890, 155.2362900},
}};

inline constexpr Medium dense_medium{ConstantIndex{1.5}};

inline constexpr std::array<MediumId, 6> medium_catalog_ids{
    MediumId::air,
    MediumId::n_bk7,
    MediumId::n_sf6,
    MediumId::n_lak9,
    MediumId::n_sf11,
    MediumId::dense,
};

[[nodiscard]] inline constexpr Medium medium(const MediumId id) {
    switch (id) {
        case MediumId::air:
            return air_medium;
        case MediumId::n_bk7:
            return n_bk7_medium;
        case MediumId::n_sf6:
            return n_sf6_medium;
        case MediumId::n_lak9:
            return n_lak9_medium;
        case MediumId::n_sf11:
            return n_sf11_medium;
        case MediumId::dense:
            return dense_medium;
    }

    return air_medium;
}

[[nodiscard]] inline constexpr std::string_view medium_name(const MediumId id) {
    switch (id) {
        case MediumId::air:
            return "Air";
        case MediumId::n_bk7:
            return "N-BK7";
        case MediumId::n_sf6:
            return "N-SF6";
        case MediumId::n_lak9:
            return "N-LAK9";
        case MediumId::n_sf11:
            return "N-SF11";
        case MediumId::dense:
            return "Dense";
    }

    return "Unknown";
}

} // namespace opsys
