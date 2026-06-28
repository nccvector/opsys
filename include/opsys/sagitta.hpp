#pragma once

#include <cmath>
#include <concepts>
#include <optional>
#include <type_traits>
#include <variant>

namespace opsys {

namespace detail {

template <std::floating_point T>
inline constexpr T k_nearly_zero = T{1.0e-14};

template <std::floating_point T>
[[nodiscard]] inline bool near_zero(T value) {
    return std::abs(value) < k_nearly_zero<T>;
}

template <std::floating_point T>
[[nodiscard]] inline std::optional<T> conic_sagitta(
    const T radius_mm,
    const T sphere_radius_mm,
    const T conic_constant) {
    if (near_zero(sphere_radius_mm)) {
        return T{0};
    }

    const T c = T{1} / sphere_radius_mm;
    const T alpha = (T{1} + conic_constant) * c * c;
    const T q2 = T{1} - alpha * radius_mm * radius_mm;
    if (q2 < T{0}) {
        return std::nullopt;
    }

    const T q = std::sqrt(q2);
    return c * radius_mm * radius_mm / (T{1} + q);
}

template <std::floating_point T>
[[nodiscard]] inline std::optional<T> conic_sagitta_derivative(
    const T radius_mm,
    const T sphere_radius_mm,
    const T conic_constant) {
    if (near_zero(sphere_radius_mm)) {
        return T{0};
    }

    const T c = T{1} / sphere_radius_mm;
    const T alpha = (T{1} + conic_constant) * c * c;
    const T q2 = T{1} - alpha * radius_mm * radius_mm;
    if (q2 <= T{0}) {
        return std::nullopt;
    }

    const T q = std::sqrt(q2);
    const T denominator = T{1} + q;
    return (T{2} * c * radius_mm * denominator + c * radius_mm * radius_mm * alpha * radius_mm / q)
        / (denominator * denominator);
}

} // namespace detail

struct PlaneSagitta {};

template <std::floating_point T>
struct ConicSagittaT {
    // Signed radius through curvature = 1 / radius. Zero curvature is a plane.
    T radius_mm{};
    T conic_constant{};
};

using ConicSagitta = ConicSagittaT<double>;

template <std::floating_point T>
struct SagittaSurfaceT {
    constexpr SagittaSurfaceT() : model(PlaneSagitta{}) {}
    explicit constexpr SagittaSurfaceT(PlaneSagitta plane) : model(plane) {}
    explicit constexpr SagittaSurfaceT(ConicSagittaT<T> conic) : model(conic) {}

    [[nodiscard]] std::optional<T> sagitta_mm(T radius_mm) const;

    [[nodiscard]] std::optional<T> dsagitta_dr(T radius_mm) const;

    std::variant<PlaneSagitta, ConicSagittaT<T>> model;
};

using SagittaSurface = SagittaSurfaceT<double>;

template <std::floating_point T>
[[nodiscard]] inline std::optional<T> sagitta_mm(const SagittaSurfaceT<T>& surface, T radius_mm) {
    return std::visit([radius_mm](const auto& model) -> std::optional<T> {
        using Model = std::decay_t<decltype(model)>;

        if constexpr (std::is_same_v<Model, PlaneSagitta>) {
            return T{0};
        } else {
            return detail::conic_sagitta(
                radius_mm,
                model.radius_mm,
                model.conic_constant);
        }
    }, surface.model);
}

template <std::floating_point T>
[[nodiscard]] inline std::optional<T> dsagitta_dr(const SagittaSurfaceT<T>& surface, T radius_mm) {
    return std::visit([radius_mm](const auto& model) -> std::optional<T> {
        using Model = std::decay_t<decltype(model)>;

        if constexpr (std::is_same_v<Model, PlaneSagitta>) {
            return T{0};
        } else {
            return detail::conic_sagitta_derivative(
                radius_mm,
                model.radius_mm,
                model.conic_constant);
        }
    }, surface.model);
}

template <std::floating_point T>
inline std::optional<T> SagittaSurfaceT<T>::sagitta_mm(T radius_mm) const {
    return opsys::sagitta_mm(*this, radius_mm);
}

template <std::floating_point T>
inline std::optional<T> SagittaSurfaceT<T>::dsagitta_dr(T radius_mm) const {
    return opsys::dsagitta_dr(*this, radius_mm);
}

} // namespace opsys
