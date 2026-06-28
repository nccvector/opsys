#pragma once

#include "osys/medium.hpp"
#include "osys/ray.hpp"
#include "osys/sagitta.hpp"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace osys {

namespace detail {

constexpr double k_ray_epsilon_mm = 1.0e-9;

[[nodiscard]] inline bool finite_ray(const Ray& ray) {
    return std::isfinite(ray.origin_mm.x)
        && std::isfinite(ray.origin_mm.y)
        && std::isfinite(ray.origin_mm.z)
        && std::isfinite(ray.direction.x)
        && std::isfinite(ray.direction.y)
        && std::isfinite(ray.direction.z)
        && std::isfinite(ray.wavelength_nm)
        && length(ray.direction) > 0.0
        && ray.wavelength_nm > 0.0;
}

[[nodiscard]] inline double radial_distance(Vec3 point) {
    return std::hypot(point.x, point.y);
}

[[nodiscard]] inline std::optional<double> forward_t(double t) {
    if (t <= k_ray_epsilon_mm || !std::isfinite(t)) {
        return std::nullopt;
    }
    return t;
}

[[nodiscard]] inline std::optional<double> smaller_forward_t(const double lhs, const double rhs) {
    const std::optional<double> forward_lhs = forward_t(lhs);
    const std::optional<double> forward_rhs = forward_t(rhs);
    if (forward_lhs.has_value() && forward_rhs.has_value()) {
        return std::min(*forward_lhs, *forward_rhs);
    }
    if (forward_lhs.has_value()) {
        return forward_lhs;
    }
    return forward_rhs;
}

[[nodiscard]] inline std::optional<double> solve_quadratic_forward(const double a, const double b, const double c) {
    if (near_zero(a)) {
        if (near_zero(b)) {
            return std::nullopt;
        }
        return forward_t(-c / b);
    }

    const double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0) {
        return std::nullopt;
    }

    const double root = std::sqrt(discriminant);
    return smaller_forward_t((-b - root) / (2.0 * a), (-b + root) / (2.0 * a));
}

[[nodiscard]] inline std::optional<Vec3> refract(
    const Vec3 &incident,
    const Vec3 &surface_normal,
    const double n_incident,
    const double n_transmitted) {
    Vec3 normal = normalize(surface_normal);
    const Vec3 direction = normalize(incident);

    if (dot(direction, normal) > 0.0) {
        normal = -normal;
    }

    const double eta = n_incident / n_transmitted;
    const double cos_i = -dot(normal, direction);
    const double k = 1.0 - eta * eta * (1.0 - cos_i * cos_i);
    if (k < 0.0) {
        return std::nullopt;
    }

    return normalize(eta * direction + (eta * cos_i - std::sqrt(k)) * normal);
}

} // namespace detail

struct OpticalSurface {
    double vertex_z_mm{};
    double aperture_radius_mm{};
    Medium medium_after{};
    SagittaSurface shape;
};

enum class TraceStatus {
    ok,
    invalid_ray,
    no_intersection,
    missed_aperture,
    total_internal_reflection,
};

struct TraceResult {
    TraceStatus status{TraceStatus::ok};
    Ray output_ray{};
    std::size_t surface_index{};
};

struct OpticalSystem {
    void add_surface(const OpticalSurface &surface);

    [[nodiscard]] TraceResult trace(const Ray& input) const;

    Medium initial_medium{};
    std::vector<OpticalSurface> surfaces;
};

inline void add_surface(OpticalSystem& system, const OpticalSurface &surface) {
    system.surfaces.push_back(surface);
}

[[nodiscard]] inline std::optional<double> intersect_surface(const Ray& ray, const OpticalSurface& surface) {
    return std::visit([&](const auto& sagitta) -> std::optional<double> {
        using Sagitta = std::decay_t<decltype(sagitta)>;

        if constexpr (std::is_same_v<Sagitta, PlaneSagitta>) {
            if (detail::near_zero(ray.direction.z)) {
                return std::nullopt;
            }
            return detail::forward_t((surface.vertex_z_mm - ray.origin_mm.z) / ray.direction.z);
        } else {
            if (detail::near_zero(sagitta.radius_mm)) {
                return std::nullopt;
            }

            const double k = 1.0 + sagitta.conic_constant;
            const double z0 = ray.origin_mm.z - surface.vertex_z_mm;
            const double a = ray.direction.x * ray.direction.x
                + ray.direction.y * ray.direction.y
                + k * ray.direction.z * ray.direction.z;
            const double b = 2.0 * (ray.origin_mm.x * ray.direction.x
                + ray.origin_mm.y * ray.direction.y
                + k * z0 * ray.direction.z
                - sagitta.radius_mm * ray.direction.z);
            const double c = ray.origin_mm.x * ray.origin_mm.x
                + ray.origin_mm.y * ray.origin_mm.y
                + k * z0 * z0
                - 2.0 * sagitta.radius_mm * z0;

            return detail::solve_quadratic_forward(a, b, c);
        }
    }, surface.shape.model);
}

[[nodiscard]] inline std::optional<Vec3> surface_normal(Vec3 point, const OpticalSurface& surface) {
    return std::visit([&](const auto& sagitta) -> std::optional<Vec3> {
        using Sagitta = std::decay_t<decltype(sagitta)>;

        if constexpr (std::is_same_v<Sagitta, PlaneSagitta>) {
            return Vec3{0.0, 0.0, 1.0};
        } else {
            if (detail::near_zero(sagitta.radius_mm)) {
                return std::nullopt;
            }

            const double z = point.z - surface.vertex_z_mm;
            const double k = 1.0 + sagitta.conic_constant;
            return normalize(Vec3{
                -point.x,
                -point.y,
                sagitta.radius_mm - k * z,
            });
        }
    }, surface.shape.model);
}

[[nodiscard]] inline TraceResult trace(const OpticalSystem& system, const Ray& input) {
    if (!detail::finite_ray(input)) {
        return {.status = TraceStatus::invalid_ray, .output_ray = input};
    }

    Ray ray = input;
    Medium current_medium = system.initial_medium;

    for (std::size_t i = 0; i < system.surfaces.size(); ++i) {
        const OpticalSurface& surface = system.surfaces[i];
        const Medium& before = current_medium;
        const Medium& after = surface.medium_after;

        const std::optional<double> hit_t = intersect_surface(ray, surface);
        if (!hit_t.has_value()) {
            return {.status = TraceStatus::no_intersection, .output_ray = ray, .surface_index = i};
        }

        const Vec3 hit = ray.origin_mm + *hit_t * ray.direction;
        if (detail::radial_distance(hit) > surface.aperture_radius_mm) {
            return {.status = TraceStatus::missed_aperture, .output_ray = ray, .surface_index = i};
        }

        const std::optional<Vec3> normal = surface_normal(hit, surface);
        if (!normal.has_value()) {
            return {.status = TraceStatus::no_intersection, .output_ray = ray, .surface_index = i};
        }

        const double n_before = refractive_index(before, ray.wavelength_nm);
        const double n_after = refractive_index(after, ray.wavelength_nm);
        const std::optional<Vec3> refracted = detail::refract(ray.direction, *normal, n_before, n_after);
        if (!refracted.has_value()) {
            return {.status = TraceStatus::total_internal_reflection, .output_ray = ray, .surface_index = i};
        }

        ray.origin_mm = hit + detail::k_ray_epsilon_mm * *refracted;
        ray.direction = *refracted;
        current_medium = surface.medium_after;
    }

    return {.status = TraceStatus::ok, .output_ray = ray, .surface_index = system.surfaces.size()};
}

inline void OpticalSystem::add_surface(const OpticalSurface &surface) {
    osys::add_surface(*this, surface);
}

inline TraceResult OpticalSystem::trace(const Ray& input) const {
    return osys::trace(*this, input);
}

} // namespace osys
