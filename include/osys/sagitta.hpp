#pragma once

#include <cmath>
#include <optional>
#include <type_traits>
#include <variant>

namespace osys {

namespace detail {

constexpr double k_nearly_zero = 1.0e-14;

[[nodiscard]] inline bool near_zero(double value) {
    return std::abs(value) < k_nearly_zero;
}

[[nodiscard]] inline std::optional<double> conic_sagitta(
    const double radius_mm,
    const double sphere_radius_mm,
    const double conic_constant) {
    if (near_zero(sphere_radius_mm)) {
        return 0.0;
    }

    const double c = 1.0 / sphere_radius_mm;
    const double alpha = (1.0 + conic_constant) * c * c;
    const double q2 = 1.0 - alpha * radius_mm * radius_mm;
    if (q2 < 0.0) {
        return std::nullopt;
    }

    const double q = std::sqrt(q2);
    return c * radius_mm * radius_mm / (1.0 + q);
}

[[nodiscard]] inline std::optional<double> conic_sagitta_derivative(
    const double radius_mm,
    const double sphere_radius_mm,
    const double conic_constant) {
    if (near_zero(sphere_radius_mm)) {
        return 0.0;
    }

    const double c = 1.0 / sphere_radius_mm;
    const double alpha = (1.0 + conic_constant) * c * c;
    const double q2 = 1.0 - alpha * radius_mm * radius_mm;
    if (q2 <= 0.0) {
        return std::nullopt;
    }

    const double q = std::sqrt(q2);
    const double denominator = 1.0 + q;
    return (2.0 * c * radius_mm * denominator + c * radius_mm * radius_mm * alpha * radius_mm / q)
        / (denominator * denominator);
}

} // namespace detail

struct PlaneSagitta {};

struct ConicSagitta {
    // Signed radius through curvature = 1 / radius. Zero curvature is a plane.
    double radius_mm{};
    double conic_constant{};
};

struct SagittaSurface {
    SagittaSurface() : model(PlaneSagitta{}) {}
    explicit SagittaSurface(PlaneSagitta plane) : model(plane) {}
    explicit SagittaSurface(ConicSagitta conic) : model(conic) {}

    [[nodiscard]] std::optional<double> sagitta_mm(double radius_mm) const;

    [[nodiscard]] std::optional<double> dsagitta_dr(double radius_mm) const;

    std::variant<PlaneSagitta, ConicSagitta> model;
};

[[nodiscard]] inline std::optional<double> sagitta_mm(const SagittaSurface& surface, double radius_mm) {
    return std::visit([radius_mm](const auto& model) -> std::optional<double> {
        using Model = std::decay_t<decltype(model)>;

        if constexpr (std::is_same_v<Model, PlaneSagitta>) {
            return 0.0;
        } else {
            return detail::conic_sagitta(
                radius_mm,
                model.radius_mm,
                model.conic_constant);
        }
    }, surface.model);
}

[[nodiscard]] inline std::optional<double> dsagitta_dr(const SagittaSurface& surface, double radius_mm) {
    return std::visit([radius_mm](const auto& model) -> std::optional<double> {
        using Model = std::decay_t<decltype(model)>;

        if constexpr (std::is_same_v<Model, PlaneSagitta>) {
            return 0.0;
        } else {
            return detail::conic_sagitta_derivative(
                radius_mm,
                model.radius_mm,
                model.conic_constant);
        }
    }, surface.model);
}

inline std::optional<double> SagittaSurface::sagitta_mm(double radius_mm) const {
    return osys::sagitta_mm(*this, radius_mm);
}

inline std::optional<double> SagittaSurface::dsagitta_dr(double radius_mm) const {
    return osys::dsagitta_dr(*this, radius_mm);
}

} // namespace osys
