#pragma once

#include <cmath>
#include <optional>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace lenses {

namespace detail {

constexpr double k_nearly_zero = 1.0e-14;

[[nodiscard]] inline bool near_zero(double value) {
    return std::abs(value) < k_nearly_zero;
}

[[nodiscard]] inline std::optional<double> spherical_sag(double radius_mm, double sphere_radius_mm) {
    if (near_zero(sphere_radius_mm)) {
        return std::nullopt;
    }

    const double under_root = sphere_radius_mm * sphere_radius_mm - radius_mm * radius_mm;
    if (under_root < 0.0) {
        return std::nullopt;
    }

    const double sign = sphere_radius_mm > 0.0 ? 1.0 : -1.0;
    return sphere_radius_mm - sign * std::sqrt(under_root);
}

[[nodiscard]] inline std::optional<double> spherical_derivative(double radius_mm, double sphere_radius_mm) {
    if (near_zero(sphere_radius_mm)) {
        return std::nullopt;
    }

    const double under_root = sphere_radius_mm * sphere_radius_mm - radius_mm * radius_mm;
    if (under_root <= 0.0) {
        return std::nullopt;
    }

    const double sign = sphere_radius_mm > 0.0 ? 1.0 : -1.0;
    return sign * radius_mm / std::sqrt(under_root);
}

[[nodiscard]] inline std::optional<double> conic_base_sag(
    double radius_mm,
    double sphere_radius_mm,
    double conic_constant) {
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

[[nodiscard]] inline std::optional<double> conic_base_derivative(
    double radius_mm,
    double sphere_radius_mm,
    double conic_constant) {
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

struct PlaneSag {};

struct SphericalSag {
    // Signed radius. R > 0 means center of curvature is in +z.
    double radius_mm{};
};

struct EvenAsphereSag {
    // Signed radius through curvature = 1 / radius. Zero curvature is a plane.
    double radius_mm{};
    double conic_constant{};
    // A4, A6, A8... coefficients in mm^(1-power).
    std::vector<double> even_coefficients;
};

struct SagSurface {
    SagSurface() : model(PlaneSag{}) {}
    explicit SagSurface(PlaneSag plane) : model(plane) {}
    explicit SagSurface(SphericalSag sphere) : model(sphere) {}
    explicit SagSurface(EvenAsphereSag asphere) : model(std::move(asphere)) {}

    [[nodiscard]] std::optional<double> sag_mm(double radius_mm) const;

    [[nodiscard]] std::optional<double> dsag_dr(double radius_mm) const;

    std::variant<PlaneSag, SphericalSag, EvenAsphereSag> model;
};

[[nodiscard]] inline std::optional<double> sag_mm(const SagSurface& surface, double radius_mm) {
    return std::visit([radius_mm](const auto& model) -> std::optional<double> {
        using Model = std::decay_t<decltype(model)>;

        if constexpr (std::is_same_v<Model, PlaneSag>) {
            return 0.0;
        } else if constexpr (std::is_same_v<Model, SphericalSag>) {
            return detail::spherical_sag(radius_mm, model.radius_mm);
        } else {
            std::optional<double> sag = detail::conic_base_sag(
                radius_mm,
                model.radius_mm,
                model.conic_constant);
            if (!sag.has_value()) {
                return std::nullopt;
            }

            double r_power = radius_mm * radius_mm;
            for (double coefficient : model.even_coefficients) {
                r_power *= radius_mm * radius_mm;
                *sag += coefficient * r_power;
            }

            return sag;
        }
    }, surface.model);
}

[[nodiscard]] inline std::optional<double> dsag_dr(const SagSurface& surface, double radius_mm) {
    return std::visit([radius_mm](const auto& model) -> std::optional<double> {
        using Model = std::decay_t<decltype(model)>;

        if constexpr (std::is_same_v<Model, PlaneSag>) {
            return 0.0;
        } else if constexpr (std::is_same_v<Model, SphericalSag>) {
            return detail::spherical_derivative(radius_mm, model.radius_mm);
        } else {
            std::optional<double> derivative = detail::conic_base_derivative(
                radius_mm,
                model.radius_mm,
                model.conic_constant);
            if (!derivative.has_value()) {
                return std::nullopt;
            }

            double r_power = radius_mm * radius_mm * radius_mm;
            double power = 4.0;
            for (double coefficient : model.even_coefficients) {
                *derivative += power * coefficient * r_power;
                r_power *= radius_mm * radius_mm;
                power += 2.0;
            }

            return derivative;
        }
    }, surface.model);
}

inline std::optional<double> SagSurface::sag_mm(double radius_mm) const {
    return lenses::sag_mm(*this, radius_mm);
}

inline std::optional<double> SagSurface::dsag_dr(double radius_mm) const {
    return lenses::dsag_dr(*this, radius_mm);
}

} // namespace lenses
