#pragma once

#include "opsys/medium.hpp"
#include "opsys/ray.hpp"
#include "opsys/sagitta.hpp"

#include <cmath>
#include <concepts>
#include <cstddef>
#include <optional>
#include <type_traits>
#include <vector>

namespace opsys {

namespace detail {

template <std::floating_point T>
inline constexpr T k_ray_epsilon_mm = T{1.0e-9};

template <std::floating_point T>
inline constexpr T k_surface_intersection_tolerance_mm = T{1.0e-7};

template <std::floating_point T>
[[nodiscard]] inline T dot(
    const T ax,
    const T ay,
    const T az,
    const T bx,
    const T by,
    const T bz) {
    return ax * bx + ay * by + az * bz;
}

template <std::floating_point T>
[[nodiscard]] inline T length(const T x, const T y, const T z) {
    return std::sqrt(dot(x, y, z, x, y, z));
}

template <std::floating_point T>
inline void normalize(T& x, T& y, T& z) {
    const T len = length(x, y, z);
    x /= len;
    y /= len;
    z /= len;
}

template <SpectralRay Ray>
[[nodiscard]] inline bool finite_ray(const Ray& ray) {
    using T = ray_scalar_t<Ray>;
    return std::isfinite(ray.ox)
        && std::isfinite(ray.oy)
        && std::isfinite(ray.oz)
        && std::isfinite(ray.dx)
        && std::isfinite(ray.dy)
        && std::isfinite(ray.dz)
        && std::isfinite(ray.wavelength)
        && length(
            static_cast<T>(ray.dx),
            static_cast<T>(ray.dy),
            static_cast<T>(ray.dz)) > T{0}
        && ray.wavelength > T{0};
}

template <std::floating_point T>
[[nodiscard]] inline T radial_distance(const T x, const T y) {
    return std::hypot(x, y);
}

template <std::floating_point T>
[[nodiscard]] inline std::optional<T> forward_t(T t) {
    if (t <= k_ray_epsilon_mm<T> || !std::isfinite(t)) {
        return std::nullopt;
    }
    return t;
}

template <std::floating_point T>
[[nodiscard]] inline bool refract(
    const T incident_x,
    const T incident_y,
    const T incident_z,
    T normal_x,
    T normal_y,
    T normal_z,
    const T n_incident,
    const T n_transmitted,
    T& refracted_x,
    T& refracted_y,
    T& refracted_z) {
    normalize(normal_x, normal_y, normal_z);

    T direction_x = incident_x;
    T direction_y = incident_y;
    T direction_z = incident_z;
    normalize(direction_x, direction_y, direction_z);

    if (dot(direction_x, direction_y, direction_z, normal_x, normal_y, normal_z) > T{0}) {
        normal_x = -normal_x;
        normal_y = -normal_y;
        normal_z = -normal_z;
    }

    const T eta = n_incident / n_transmitted;
    const T cos_i = -dot(normal_x, normal_y, normal_z, direction_x, direction_y, direction_z);
    const T k = T{1} - eta * eta * (T{1} - cos_i * cos_i);
    if (k < T{0}) {
        return false;
    }

    refracted_x = eta * direction_x + (eta * cos_i - std::sqrt(k)) * normal_x;
    refracted_y = eta * direction_y + (eta * cos_i - std::sqrt(k)) * normal_y;
    refracted_z = eta * direction_z + (eta * cos_i - std::sqrt(k)) * normal_z;
    normalize(refracted_x, refracted_y, refracted_z);
    return true;
}

} // namespace detail

template <std::floating_point T>
struct OpticalSurfaceT {
    T vertex_z_mm{};
    T aperture_radius_mm{};
    MediumT<T> medium_after{};
    SagittaSurfaceT<T> shape;
};

using OpticalSurface = OpticalSurfaceT<double>;

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

template <std::floating_point T>
struct OpticalSystemT {
    void add_surface(const OpticalSurfaceT<T>& surface);

    template <SpectralRay Ray>
    [[nodiscard]] TraceResult<Ray> trace(const Ray& input) const;

    template <SpectralRay Ray>
    [[nodiscard]] TraceResult<Ray> reverse_trace(const Ray& input) const;

    MediumT<T> initial_medium{};
    std::vector<OpticalSurfaceT<T>> surfaces;
};

using OpticalSystem = OpticalSystemT<double>;

namespace detail {

template <std::floating_point T>
[[nodiscard]] inline const MediumT<T>& medium_before_surface(const OpticalSystemT<T>& system, const std::size_t surface_index) {
    if (surface_index == 0) {
        return system.initial_medium;
    }

    return system.surfaces[surface_index - 1].medium_after;
}

template <std::floating_point T>
[[nodiscard]] inline const MediumT<T>& image_side_medium(const OpticalSystemT<T>& system) {
    if (system.surfaces.empty()) {
        return system.initial_medium;
    }

    return system.surfaces.back().medium_after;
}

} // namespace detail

template <std::floating_point T>
inline void add_surface(OpticalSystemT<T>& system, const OpticalSurfaceT<T>& surface) {
    system.surfaces.push_back(surface);
}

template <std::floating_point T, SpectralRay Ray>
[[nodiscard]] inline std::optional<T> intersect_surface(const Ray& ray, const OpticalSurfaceT<T>& surface) {
    return std::visit([&](const auto& sagitta) -> std::optional<T> {
        using Sagitta = std::decay_t<decltype(sagitta)>;

        if constexpr (std::is_same_v<Sagitta, PlaneSagitta>) {
            const T ray_dz = static_cast<T>(ray.dz);
            if (detail::near_zero(ray_dz)) {
                return std::nullopt;
            }
            return detail::forward_t((surface.vertex_z_mm - static_cast<T>(ray.oz)) / ray_dz);
        } else {
            if (detail::near_zero(sagitta.radius_mm)) {
                return std::nullopt;
            }

            const T ray_ox = static_cast<T>(ray.ox);
            const T ray_oy = static_cast<T>(ray.oy);
            const T ray_oz = static_cast<T>(ray.oz);
            const T ray_dx = static_cast<T>(ray.dx);
            const T ray_dy = static_cast<T>(ray.dy);
            const T ray_dz = static_cast<T>(ray.dz);
            const T k = T{1} + sagitta.conic_constant;
            const T z0 = ray_oz - surface.vertex_z_mm;
            const T a = ray_dx * ray_dx + ray_dy * ray_dy + k * ray_dz * ray_dz;
            const T b = T{2} * (ray_ox * ray_dx
                + ray_oy * ray_dy
                + k * z0 * ray_dz
                - sagitta.radius_mm * ray_dz);
            const T c = ray_ox * ray_ox + ray_oy * ray_oy + k * z0 * z0 - T{2} * sagitta.radius_mm * z0;

            std::optional<T> best_t;
            const auto consider_root = [&](const T t) {
                const std::optional<T> forward = detail::forward_t(t);
                if (!forward.has_value()) {
                    return;
                }

                const T hit_x = ray_ox + *forward * ray_dx;
                const T hit_y = ray_oy + *forward * ray_dy;
                const T hit_z = ray_oz + *forward * ray_dz;
                const std::optional<T> expected_z = detail::conic_sagitta(
                    detail::radial_distance(hit_x, hit_y),
                    sagitta.radius_mm,
                    sagitta.conic_constant);
                if (!expected_z.has_value()) {
                    return;
                }

                const T local_z = hit_z - surface.vertex_z_mm;
                if (std::abs(local_z - *expected_z) > detail::k_surface_intersection_tolerance_mm<T>) {
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
                const T discriminant = b * b - T{4} * a * c;
                if (discriminant < T{0}) {
                    return std::nullopt;
                }

                const T root = std::sqrt(discriminant);
                consider_root((-b - root) / (T{2} * a));
                consider_root((-b + root) / (T{2} * a));
            }

            return best_t;
        }
    }, surface.shape.model);
}

template <std::floating_point T>
[[nodiscard]] inline bool surface_normal(
    const T point_x,
    const T point_y,
    const T point_z,
    const OpticalSurfaceT<T>& surface,
    T& normal_x,
    T& normal_y,
    T& normal_z) {
    return std::visit([&](const auto& sagitta) -> bool {
        using Sagitta = std::decay_t<decltype(sagitta)>;

        if constexpr (std::is_same_v<Sagitta, PlaneSagitta>) {
            normal_x = T{0};
            normal_y = T{0};
            normal_z = T{1};
            return true;
        } else {
            if (detail::near_zero(sagitta.radius_mm)) {
                return false;
            }

            const T z = point_z - surface.vertex_z_mm;
            const T k = T{1} + sagitta.conic_constant;
            normal_x = -point_x;
            normal_y = -point_y;
            normal_z = sagitta.radius_mm - k * z;
            detail::normalize(normal_x, normal_y, normal_z);
            return true;
        }
    }, surface.shape.model);
}

namespace detail {

template <std::floating_point T, SpectralRay Ray>
[[nodiscard]] inline TraceStatus trace_surface(
    Ray& ray,
    MediumT<T>& current_medium,
    const OpticalSurfaceT<T>& surface,
    const MediumT<T>& next_medium) {
    const std::optional<T> hit_t = intersect_surface(ray, surface);
    if (!hit_t.has_value()) {
        return TraceStatus::no_intersection;
    }

    const T hit_x = static_cast<T>(ray.ox) + *hit_t * static_cast<T>(ray.dx);
    const T hit_y = static_cast<T>(ray.oy) + *hit_t * static_cast<T>(ray.dy);
    const T hit_z = static_cast<T>(ray.oz) + *hit_t * static_cast<T>(ray.dz);
    if (radial_distance(hit_x, hit_y) > surface.aperture_radius_mm) {
        return TraceStatus::missed_aperture;
    }

    T normal_x{};
    T normal_y{};
    T normal_z{};
    if (!surface_normal(hit_x, hit_y, hit_z, surface, normal_x, normal_y, normal_z)) {
        return TraceStatus::no_intersection;
    }

    const T wavelength = static_cast<T>(ray.wavelength);
    const T n_before = refractive_index(current_medium, wavelength);
    const T n_after = refractive_index(next_medium, wavelength);
    T refracted_x{};
    T refracted_y{};
    T refracted_z{};
    if (!refract(
            static_cast<T>(ray.dx),
            static_cast<T>(ray.dy),
            static_cast<T>(ray.dz),
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

    ray.ox = hit_x + k_ray_epsilon_mm<T> * refracted_x;
    ray.oy = hit_y + k_ray_epsilon_mm<T> * refracted_y;
    ray.oz = hit_z + k_ray_epsilon_mm<T> * refracted_z;
    ray.dx = refracted_x;
    ray.dy = refracted_y;
    ray.dz = refracted_z;
    current_medium = next_medium;
    return TraceStatus::ok;
}

} // namespace detail

template <std::floating_point T, SpectralRay Ray>
[[nodiscard]] inline TraceResult<Ray> trace(const OpticalSystemT<T>& system, const Ray& input) {
    if (!detail::finite_ray(input)) {
        return {.status = TraceStatus::invalid_ray, .output_ray = input};
    }

    Ray ray = input;
    MediumT<T> current_medium = system.initial_medium;

    for (std::size_t i = 0; i < system.surfaces.size(); ++i) {
        const OpticalSurfaceT<T>& surface = system.surfaces[i];
        const TraceStatus status = detail::trace_surface(ray, current_medium, surface, surface.medium_after);
        if (status != TraceStatus::ok) {
            return {.status = status, .output_ray = ray, .surface_index = i};
        }
    }

    return {.status = TraceStatus::ok, .output_ray = ray, .surface_index = system.surfaces.size()};
}

template <std::floating_point T, SpectralRay Ray>
[[nodiscard]] inline TraceResult<Ray> reverse_trace(const OpticalSystemT<T>& system, const Ray& input) {
    if (!detail::finite_ray(input)) {
        return {.status = TraceStatus::invalid_ray, .output_ray = input};
    }

    Ray ray = input;
    MediumT<T> current_medium = detail::image_side_medium(system);

    for (std::size_t remaining = system.surfaces.size(); remaining > 0; --remaining) {
        const std::size_t i = remaining - 1;
        const OpticalSurfaceT<T>& surface = system.surfaces[i];
        const MediumT<T>& next_medium = detail::medium_before_surface(system, i);
        const TraceStatus status = detail::trace_surface(ray, current_medium, surface, next_medium);
        if (status != TraceStatus::ok) {
            return {.status = status, .output_ray = ray, .surface_index = i};
        }
    }

    return {.status = TraceStatus::ok, .output_ray = ray, .surface_index = system.surfaces.size()};
}

template <std::floating_point T>
inline void OpticalSystemT<T>::add_surface(const OpticalSurfaceT<T>& surface) {
    opsys::add_surface(*this, surface);
}

template <std::floating_point T>
template <SpectralRay Ray>
inline TraceResult<Ray> OpticalSystemT<T>::trace(const Ray& input) const {
    return opsys::trace(*this, input);
}

template <std::floating_point T>
template <SpectralRay Ray>
inline TraceResult<Ray> OpticalSystemT<T>::reverse_trace(const Ray& input) const {
    return opsys::reverse_trace(*this, input);
}

} // namespace opsys
