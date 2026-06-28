#pragma once

#include "opsys/medium.hpp"
#include "opsys/ray.hpp"
#include "opsys/sagitta.hpp"

#include <cstddef>
#include <cmath>
#include <optional>
#include <type_traits>
#include <vector>

namespace opsys {

namespace detail {

constexpr double k_ray_epsilon_mm = 1.0e-9;
constexpr double k_surface_intersection_tolerance_mm = 1.0e-7;

[[nodiscard]] inline bool finite_ray(const RayLike auto& ray) {
    const Vec3 origin = ray_origin_mm(ray);
    const Vec3 direction = ray_direction(ray);
    const double wavelength = ray_wavelength(ray);
    return std::isfinite(origin.x)
        && std::isfinite(origin.y)
        && std::isfinite(origin.z)
        && std::isfinite(direction.x)
        && std::isfinite(direction.y)
        && std::isfinite(direction.z)
        && std::isfinite(wavelength)
        && length(direction) > 0.0
        && wavelength > 0.0;
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

template <RayLike Ray>
struct TraceResult {
    TraceStatus status{TraceStatus::ok};
    Ray output_ray;
    std::size_t surface_index{};
};

struct OpticalSystem {
    void add_surface(const OpticalSurface &surface);

    template <RayLike Ray>
    [[nodiscard]] TraceResult<Ray> trace(const Ray& input) const;

    template <RayLike Ray>
    [[nodiscard]] TraceResult<Ray> reverse_trace(const Ray& input) const;

    Medium initial_medium{};
    std::vector<OpticalSurface> surfaces;
};

namespace detail {

[[nodiscard]] inline const Medium& medium_before_surface(const OpticalSystem& system, const std::size_t surface_index) {
    if (surface_index == 0) {
        return system.initial_medium;
    }

    return system.surfaces[surface_index - 1].medium_after;
}

[[nodiscard]] inline const Medium& image_side_medium(const OpticalSystem& system) {
    if (system.surfaces.empty()) {
        return system.initial_medium;
    }

    return system.surfaces.back().medium_after;
}

} // namespace detail

inline void add_surface(OpticalSystem& system, const OpticalSurface &surface) {
    system.surfaces.push_back(surface);
}

[[nodiscard]] inline std::optional<double> intersect_surface(const RayLike auto& ray, const OpticalSurface& surface) {
    const Vec3 origin = ray_origin_mm(ray);
    const Vec3 direction = ray_direction(ray);
    return std::visit([&](const auto& sagitta) -> std::optional<double> {
        using Sagitta = std::decay_t<decltype(sagitta)>;

        if constexpr (std::is_same_v<Sagitta, PlaneSagitta>) {
            if (detail::near_zero(direction.z)) {
                return std::nullopt;
            }
            return detail::forward_t((surface.vertex_z_mm - origin.z) / direction.z);
        } else {
            if (detail::near_zero(sagitta.radius_mm)) {
                return std::nullopt;
            }

            const double k = 1.0 + sagitta.conic_constant;
            const double z0 = origin.z - surface.vertex_z_mm;
            const double a = direction.x * direction.x
                + direction.y * direction.y
                + k * direction.z * direction.z;
            const double b = 2.0 * (origin.x * direction.x
                + origin.y * direction.y
                + k * z0 * direction.z
                - sagitta.radius_mm * direction.z);
            const double c = origin.x * origin.x
                + origin.y * origin.y
                + k * z0 * z0
                - 2.0 * sagitta.radius_mm * z0;

            std::optional<double> best_t;
            const auto consider_root = [&](const double t) {
                const std::optional<double> forward = detail::forward_t(t);
                if (!forward.has_value()) {
                    return;
                }

                const Vec3 point = origin + *forward * direction;
                const std::optional<double> expected_z = detail::conic_sagitta(
                    detail::radial_distance(point),
                    sagitta.radius_mm,
                    sagitta.conic_constant);
                if (!expected_z.has_value()) {
                    return;
                }

                const double local_z = point.z - surface.vertex_z_mm;
                if (std::abs(local_z - *expected_z) > detail::k_surface_intersection_tolerance_mm) {
                    return;
                }

                if (!best_t.has_value() || *forward < *best_t) {
                    best_t = forward;
                }
            };

            if (detail::near_zero(a)) {
                if (detail::near_zero(b)) {
                    return std::nullopt;
                }
                consider_root(-c / b);
            } else {
                const double discriminant = b * b - 4.0 * a * c;
                if (discriminant < 0.0) {
                    return std::nullopt;
                }

                const double root = std::sqrt(discriminant);
                consider_root((-b - root) / (2.0 * a));
                consider_root((-b + root) / (2.0 * a));
            }

            return best_t;
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

namespace detail {

[[nodiscard]] inline TraceStatus trace_surface(
    RayLike auto& ray,
    Medium& current_medium,
    const OpticalSurface& surface,
    const Medium& next_medium) {
    const std::optional<double> hit_t = intersect_surface(ray, surface);
    if (!hit_t.has_value()) {
        return TraceStatus::no_intersection;
    }

    const Vec3 origin = ray_origin_mm(ray);
    const Vec3 direction = ray_direction(ray);
    const Vec3 hit = origin + *hit_t * direction;
    if (radial_distance(hit) > surface.aperture_radius_mm) {
        return TraceStatus::missed_aperture;
    }

    const std::optional<Vec3> normal = surface_normal(hit, surface);
    if (!normal.has_value()) {
        return TraceStatus::no_intersection;
    }

    const double wavelength = ray_wavelength(ray);
    const double n_before = refractive_index(current_medium, wavelength);
    const double n_after = refractive_index(next_medium, wavelength);
    const std::optional<Vec3> refracted = refract(direction, *normal, n_before, n_after);
    if (!refracted.has_value()) {
        return TraceStatus::total_internal_reflection;
    }

    set_ray_origin_mm(ray, hit + k_ray_epsilon_mm * *refracted);
    set_ray_direction(ray, *refracted);
    current_medium = next_medium;
    return TraceStatus::ok;
}

} // namespace detail

template <RayLike Ray>
[[nodiscard]] inline TraceResult<Ray> trace(const OpticalSystem& system, const Ray& input) {
    if (!detail::finite_ray(input)) {
        return {.status = TraceStatus::invalid_ray, .output_ray = input};
    }

    Ray ray = input;
    Medium current_medium = system.initial_medium;

    for (std::size_t i = 0; i < system.surfaces.size(); ++i) {
        const OpticalSurface& surface = system.surfaces[i];
        const TraceStatus status = detail::trace_surface(ray, current_medium, surface, surface.medium_after);
        if (status != TraceStatus::ok) {
            return {.status = status, .output_ray = ray, .surface_index = i};
        }
    }

    return {.status = TraceStatus::ok, .output_ray = ray, .surface_index = system.surfaces.size()};
}

template <RayLike Ray>
[[nodiscard]] inline TraceResult<Ray> reverse_trace(const OpticalSystem& system, const Ray& input) {
    if (!detail::finite_ray(input)) {
        return {.status = TraceStatus::invalid_ray, .output_ray = input};
    }

    Ray ray = input;
    Medium current_medium = detail::image_side_medium(system);

    for (std::size_t remaining = system.surfaces.size(); remaining > 0; --remaining) {
        const std::size_t i = remaining - 1;
        const OpticalSurface& surface = system.surfaces[i];
        const Medium& next_medium = detail::medium_before_surface(system, i);
        const TraceStatus status = detail::trace_surface(ray, current_medium, surface, next_medium);
        if (status != TraceStatus::ok) {
            return {.status = status, .output_ray = ray, .surface_index = i};
        }
    }

    return {.status = TraceStatus::ok, .output_ray = ray, .surface_index = system.surfaces.size()};
}

inline void OpticalSystem::add_surface(const OpticalSurface &surface) {
    opsys::add_surface(*this, surface);
}

template <RayLike Ray>
inline TraceResult<Ray> OpticalSystem::trace(const Ray& input) const {
    return opsys::trace(*this, input);
}

template <RayLike Ray>
inline TraceResult<Ray> OpticalSystem::reverse_trace(const Ray& input) const {
    return opsys::reverse_trace(*this, input);
}

} // namespace opsys
