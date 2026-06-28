#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <string_view>
#include <type_traits>
#include <variant>

namespace opsys {

template <std::floating_point T>
struct ConstantIndexT {
    T n{};
};

template <std::floating_point T>
struct Sellmeier3T {
    // C coefficients are in micrometer^2. Input wavelengths are still nm.
    std::array<T, 3> b{};
    std::array<T, 3> c_um2{};
};

template <std::floating_point T>
struct AbbeVdT {
    T nd{};
    T vd{};
};

using ConstantIndex = ConstantIndexT<double>;
using Sellmeier3 = Sellmeier3T<double>;
using AbbeVd = AbbeVdT<double>;

template <std::floating_point T>
struct MediumT {
    constexpr MediumT() = default;
    explicit constexpr MediumT(ConstantIndexT<T> constant) : model(constant) {}
    explicit constexpr MediumT(Sellmeier3T<T> sellmeier) : model(sellmeier) {}
    explicit constexpr MediumT(AbbeVdT<T> abbe) : model(abbe) {}

    [[nodiscard]] T refractive_index(T wavelength_nm) const;

    std::variant<ConstantIndexT<T>, Sellmeier3T<T>, AbbeVdT<T>> model{ConstantIndexT<T>{T{1}}};
};

using Medium = MediumT<double>;

template <std::floating_point T>
[[nodiscard]] inline T refractive_index(const MediumT<T>& medium, T wavelength_nm) {
    return std::visit([wavelength_nm](const auto& model) {
        using Model = std::decay_t<decltype(model)>;

        if constexpr (std::is_same_v<Model, ConstantIndexT<T>>) {
            return model.n;
        } else if constexpr (std::is_same_v<Model, Sellmeier3T<T>>) {
            const T lambda_um = wavelength_nm * T{0.001};
            const T lambda2 = lambda_um * lambda_um;
            T n2 = T{1};

            for (std::size_t i = 0; i < model.b.size(); ++i) {
                n2 += model.b[i] * lambda2 / (lambda2 - model.c_um2[i]);
            }

            return std::sqrt(n2);
        } else {
            if (model.vd <= T{0}) {
                return model.nd;
            }

            constexpr T lambda_c_um = T{0.6562725};
            constexpr T lambda_d_um = T{0.5875618};
            constexpr T lambda_f_um = T{0.4861327};
            const T lambda_um = wavelength_nm * T{0.001};
            const T delta_fc = (model.nd - T{1}) / model.vd;
            const T cauchy_b = delta_fc
                / (T{1} / (lambda_f_um * lambda_f_um) - T{1} / (lambda_c_um * lambda_c_um));
            const T cauchy_a = model.nd - cauchy_b / (lambda_d_um * lambda_d_um);

            return cauchy_a + cauchy_b / (lambda_um * lambda_um);
        }
    }, medium.model);
}

template <std::floating_point T>
inline T MediumT<T>::refractive_index(T wavelength_nm) const {
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

template <std::floating_point T>
inline constexpr MediumT<T> air_medium_v{ConstantIndexT<T>{T{1}}};

template <std::floating_point T>
inline constexpr MediumT<T> n_bk7_medium_v{Sellmeier3T<T>{
        .b = {T{1.03961212}, T{0.231792344}, T{1.01046945}},
        .c_um2 = {T{0.00600069867}, T{0.0200179144}, T{103.560653}},
}};

template <std::floating_point T>
inline constexpr MediumT<T> n_sf6_medium_v{Sellmeier3T<T>{
        .b = {T{1.72448482}, T{0.390104889}, T{1.04572858}},
        .c_um2 = {T{0.01349871764}, T{0.05694099484}, T{118.5571585}},
}};

template <std::floating_point T>
inline constexpr MediumT<T> n_lak9_medium_v{Sellmeier3T<T>{
        .b = {T{1.46231905}, T{0.344399589}, T{1.15508372}},
        .c_um2 = {T{0.00724270156}, T{0.0243353131}, T{85.4686868}},
}};

template <std::floating_point T>
inline constexpr MediumT<T> n_sf11_medium_v{Sellmeier3T<T>{
        .b = {T{1.73759695}, T{0.313747346}, T{1.89878101}},
        .c_um2 = {T{0.01318871330}, T{0.06238975890}, T{155.2362900}},
}};

template <std::floating_point T>
inline constexpr MediumT<T> dense_medium_v{ConstantIndexT<T>{T{1.5}}};

inline constexpr Medium air_medium = air_medium_v<double>;
inline constexpr Medium n_bk7_medium = n_bk7_medium_v<double>;
inline constexpr Medium n_sf6_medium = n_sf6_medium_v<double>;
inline constexpr Medium n_lak9_medium = n_lak9_medium_v<double>;
inline constexpr Medium n_sf11_medium = n_sf11_medium_v<double>;
inline constexpr Medium dense_medium = dense_medium_v<double>;

inline constexpr std::array<MediumId, 6> medium_catalog_ids{
    MediumId::air,
    MediumId::n_bk7,
    MediumId::n_sf6,
    MediumId::n_lak9,
    MediumId::n_sf11,
    MediumId::dense,
};

template <std::floating_point T = double>
[[nodiscard]] inline constexpr MediumT<T> medium(const MediumId id) {
    switch (id) {
        case MediumId::air:
            return air_medium_v<T>;
        case MediumId::n_bk7:
            return n_bk7_medium_v<T>;
        case MediumId::n_sf6:
            return n_sf6_medium_v<T>;
        case MediumId::n_lak9:
            return n_lak9_medium_v<T>;
        case MediumId::n_sf11:
            return n_sf11_medium_v<T>;
        case MediumId::dense:
            return dense_medium_v<T>;
    }

    return air_medium_v<T>;
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
