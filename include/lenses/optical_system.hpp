#pragma once

#include "lenses/material.hpp"
#include "lenses/ray.hpp"
#include "lenses/sag_surface.hpp"

#include <cstddef>
#include <cmath>
#include <optional>
#include <utility>
#include <vector>

namespace lenses {

namespace detail {

constexpr double k_ray_epsilon_mm = 1.0e-9;
constexpr double k_root_tolerance_mm = 1.0e-10;
constexpr int k_newton_iterations = 32;

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

[[nodiscard]] inline std::optional<Vec3> refract(
    Vec3 incident,
    Vec3 surface_normal,
    double n_incident,
    double n_transmitted) {
    Vec3 normal = normalize(surface_normal);
    Vec3 direction = normalize(incident);

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
    SagSurface shape;
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
    void add_surface(OpticalSurface surface);

    [[nodiscard]] TraceResult trace(const Ray& input) const;

    Medium initial_medium{};
    std::vector<OpticalSurface> surfaces;
};

inline void add_surface(OpticalSystem& system, OpticalSurface surface) {
    system.surfaces.push_back(std::move(surface));
}

[[nodiscard]] inline std::optional<double> intersect_surface(const Ray& ray, const OpticalSurface& surface) {
    if (std::abs(ray.direction.z) < 1.0e-14) {
        return std::nullopt;
    }

    double t = (surface.vertex_z_mm - ray.origin_mm.z) / ray.direction.z;
    if (t <= detail::k_ray_epsilon_mm) {
        t = detail::k_ray_epsilon_mm;
    }

    for (int iteration = 0; iteration < detail::k_newton_iterations; ++iteration) {
        const Vec3 point = ray.origin_mm + t * ray.direction;
        const double r = detail::radial_distance(point);
        const std::optional<double> sag = sag_mm(surface.shape, r);
        const std::optional<double> slope = dsag_dr(surface.shape, r);
        if (!sag.has_value() || !slope.has_value()) {
            return std::nullopt;
        }

        const double f = point.z - surface.vertex_z_mm - *sag;
        if (std::abs(f) < detail::k_root_tolerance_mm) {
            return t > detail::k_ray_epsilon_mm ? std::optional<double>{t} : std::nullopt;
        }

        double dr_dt = 0.0;
        if (r > 0.0) {
            dr_dt = (point.x * ray.direction.x + point.y * ray.direction.y) / r;
        }

        const double df_dt = ray.direction.z - *slope * dr_dt;
        if (std::abs(df_dt) < 1.0e-14) {
            return std::nullopt;
        }

        t -= f / df_dt;
        if (t <= detail::k_ray_epsilon_mm || !std::isfinite(t)) {
            return std::nullopt;
        }
    }

    return std::nullopt;
}

[[nodiscard]] inline std::optional<Vec3> surface_normal(Vec3 point, const OpticalSurface& surface) {
    const double r = detail::radial_distance(point);
    const std::optional<double> slope = dsag_dr(surface.shape, r);
    if (!slope.has_value()) {
        return std::nullopt;
    }

    if (r == 0.0) {
        return Vec3{0.0, 0.0, 1.0};
    }

    return normalize(Vec3{
        -*slope * point.x / r,
        -*slope * point.y / r,
        1.0,
    });
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

inline void OpticalSystem::add_surface(OpticalSurface surface) {
    lenses::add_surface(*this, std::move(surface));
}

inline TraceResult OpticalSystem::trace(const Ray& input) const {
    return lenses::trace(*this, input);
}

} // namespace lenses
