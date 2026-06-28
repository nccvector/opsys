#pragma once

#include "opsys/medium.hpp"
#include "opsys/ray.hpp"
#include "opsys/sagitta.hpp"

#include <cmath>
#include <cstddef>
#include <optional>
#include <type_traits>
#include <vector>

namespace opsys {

namespace detail {

constexpr double k_ray_epsilon_mm = 1.0e-9;
constexpr double k_surface_intersection_tolerance_mm = 1.0e-7;

[[nodiscard]] inline double dot(
    const double ax,
    const double ay,
    const double az,
    const double bx,
    const double by,
    const double bz) {
    return ax * bx + ay * by + az * bz;
}

[[nodiscard]] inline double length(const double x, const double y, const double z) {
    return std::sqrt(dot(x, y, z, x, y, z));
}

inline void normalize(double& x, double& y, double& z) {
    const double len = length(x, y, z);
    x /= len;
    y /= len;
    z /= len;
}

[[nodiscard]] inline bool finite_ray(const SpectralRay auto& ray) {
    return std::isfinite(ray.ox)
        && std::isfinite(ray.oy)
        && std::isfinite(ray.oz)
        && std::isfinite(ray.dx)
        && std::isfinite(ray.dy)
        && std::isfinite(ray.dz)
        && std::isfinite(ray.wavelength)
        && length(ray.dx, ray.dy, ray.dz) > 0.0
        && ray.wavelength > 0.0;
}

[[nodiscard]] inline double radial_distance(const double x, const double y) {
    return std::hypot(x, y);
}

[[nodiscard]] inline std::optional<double> forward_t(double t) {
    if (t <= k_ray_epsilon_mm || !std::isfinite(t)) {
        return std::nullopt;
    }
    return t;
}

[[nodiscard]] inline bool refract(
    const double incident_x,
    const double incident_y,
    const double incident_z,
    double normal_x,
    double normal_y,
    double normal_z,
    const double n_incident,
    const double n_transmitted,
    double& refracted_x,
    double& refracted_y,
    double& refracted_z) {
    normalize(normal_x, normal_y, normal_z);

    double direction_x = incident_x;
    double direction_y = incident_y;
    double direction_z = incident_z;
    normalize(direction_x, direction_y, direction_z);

    if (dot(direction_x, direction_y, direction_z, normal_x, normal_y, normal_z) > 0.0) {
        normal_x = -normal_x;
        normal_y = -normal_y;
        normal_z = -normal_z;
    }

    const double eta = n_incident / n_transmitted;
    const double cos_i = -dot(normal_x, normal_y, normal_z, direction_x, direction_y, direction_z);
    const double k = 1.0 - eta * eta * (1.0 - cos_i * cos_i);
    if (k < 0.0) {
        return false;
    }

    refracted_x = eta * direction_x + (eta * cos_i - std::sqrt(k)) * normal_x;
    refracted_y = eta * direction_y + (eta * cos_i - std::sqrt(k)) * normal_y;
    refracted_z = eta * direction_z + (eta * cos_i - std::sqrt(k)) * normal_z;
    normalize(refracted_x, refracted_y, refracted_z);
    return true;
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

template <SpectralRay Ray>
struct TraceResult {
    TraceStatus status{TraceStatus::ok};
    Ray output_ray;
    std::size_t surface_index{};
};

struct OpticalSystem {
    void add_surface(const OpticalSurface &surface);

    template <SpectralRay Ray>
    [[nodiscard]] TraceResult<Ray> trace(const Ray& input) const;

    template <SpectralRay Ray>
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

[[nodiscard]] inline std::optional<double> intersect_surface(const SpectralRay auto& ray, const OpticalSurface& surface) {
    return std::visit([&](const auto& sagitta) -> std::optional<double> {
        using Sagitta = std::decay_t<decltype(sagitta)>;

        if constexpr (std::is_same_v<Sagitta, PlaneSagitta>) {
            if (detail::near_zero(ray.dz)) {
                return std::nullopt;
            }
            return detail::forward_t((surface.vertex_z_mm - ray.oz) / ray.dz);
        } else {
            if (detail::near_zero(sagitta.radius_mm)) {
                return std::nullopt;
            }

            const double k = 1.0 + sagitta.conic_constant;
            const double z0 = ray.oz - surface.vertex_z_mm;
            const double a = ray.dx * ray.dx + ray.dy * ray.dy + k * ray.dz * ray.dz;
            const double b = 2.0 * (ray.ox * ray.dx
                + ray.oy * ray.dy
                + k * z0 * ray.dz
                - sagitta.radius_mm * ray.dz);
            const double c = ray.ox * ray.ox + ray.oy * ray.oy + k * z0 * z0 - 2.0 * sagitta.radius_mm * z0;

            std::optional<double> best_t;
            const auto consider_root = [&](const double t) {
                const std::optional<double> forward = detail::forward_t(t);
                if (!forward.has_value()) {
                    return;
                }

                const double hit_x = ray.ox + *forward * ray.dx;
                const double hit_y = ray.oy + *forward * ray.dy;
                const double hit_z = ray.oz + *forward * ray.dz;
                const std::optional<double> expected_z = detail::conic_sagitta(
                    detail::radial_distance(hit_x, hit_y),
                    sagitta.radius_mm,
                    sagitta.conic_constant);
                if (!expected_z.has_value()) {
                    return;
                }

                const double local_z = hit_z - surface.vertex_z_mm;
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

[[nodiscard]] inline bool surface_normal(
    const double point_x,
    const double point_y,
    const double point_z,
    const OpticalSurface& surface,
    double& normal_x,
    double& normal_y,
    double& normal_z) {
    return std::visit([&](const auto& sagitta) -> bool {
        using Sagitta = std::decay_t<decltype(sagitta)>;

        if constexpr (std::is_same_v<Sagitta, PlaneSagitta>) {
            normal_x = 0.0;
            normal_y = 0.0;
            normal_z = 1.0;
            return true;
        } else {
            if (detail::near_zero(sagitta.radius_mm)) {
                return false;
            }

            const double z = point_z - surface.vertex_z_mm;
            const double k = 1.0 + sagitta.conic_constant;
            normal_x = -point_x;
            normal_y = -point_y;
            normal_z = sagitta.radius_mm - k * z;
            detail::normalize(normal_x, normal_y, normal_z);
            return true;
        }
    }, surface.shape.model);
}

namespace detail {

[[nodiscard]] inline TraceStatus trace_surface(
    SpectralRay auto& ray,
    Medium& current_medium,
    const OpticalSurface& surface,
    const Medium& next_medium) {
    const std::optional<double> hit_t = intersect_surface(ray, surface);
    if (!hit_t.has_value()) {
        return TraceStatus::no_intersection;
    }

    const double hit_x = ray.ox + *hit_t * ray.dx;
    const double hit_y = ray.oy + *hit_t * ray.dy;
    const double hit_z = ray.oz + *hit_t * ray.dz;
    if (radial_distance(hit_x, hit_y) > surface.aperture_radius_mm) {
        return TraceStatus::missed_aperture;
    }

    double normal_x{};
    double normal_y{};
    double normal_z{};
    if (!surface_normal(hit_x, hit_y, hit_z, surface, normal_x, normal_y, normal_z)) {
        return TraceStatus::no_intersection;
    }

    const double n_before = refractive_index(current_medium, ray.wavelength);
    const double n_after = refractive_index(next_medium, ray.wavelength);
    double refracted_x{};
    double refracted_y{};
    double refracted_z{};
    if (!refract(
            ray.dx,
            ray.dy,
            ray.dz,
            normal_x,
            normal_y,
            normal_z,
            n_before,
            n_after,
            refracted_x,
            refracted_y,
            refracted_z)) {
        return TraceStatus::total_internal_reflection;
    }

    ray.ox = hit_x + k_ray_epsilon_mm * refracted_x;
    ray.oy = hit_y + k_ray_epsilon_mm * refracted_y;
    ray.oz = hit_z + k_ray_epsilon_mm * refracted_z;
    ray.dx = refracted_x;
    ray.dy = refracted_y;
    ray.dz = refracted_z;
    current_medium = next_medium;
    return TraceStatus::ok;
}

} // namespace detail

template <SpectralRay Ray>
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

template <SpectralRay Ray>
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

template <SpectralRay Ray>
inline TraceResult<Ray> OpticalSystem::trace(const Ray& input) const {
    return opsys::trace(*this, input);
}

template <SpectralRay Ray>
inline TraceResult<Ray> OpticalSystem::reverse_trace(const Ray& input) const {
    return opsys::reverse_trace(*this, input);
}

} // namespace opsys
