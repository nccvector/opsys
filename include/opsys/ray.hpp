#pragma once

#include "opsys/geometry.hpp"

#include <concepts>

namespace opsys {

template <typename Ray>
concept RayLike = requires(Ray ray, const Ray const_ray, double value) {
    { const_ray.ox } -> std::convertible_to<double>;
    { const_ray.oy } -> std::convertible_to<double>;
    { const_ray.oz } -> std::convertible_to<double>;
    { const_ray.dx } -> std::convertible_to<double>;
    { const_ray.dy } -> std::convertible_to<double>;
    { const_ray.dz } -> std::convertible_to<double>;
    { const_ray.wavelength } -> std::convertible_to<double>;

    ray.ox = value;
    ray.oy = value;
    ray.oz = value;
    ray.dx = value;
    ray.dy = value;
    ray.dz = value;
    ray.wavelength = value;
};

[[nodiscard]] constexpr Vec3 ray_origin_mm(const RayLike auto& ray) {
    return Vec3{
        static_cast<double>(ray.ox),
        static_cast<double>(ray.oy),
        static_cast<double>(ray.oz),
    };
}

[[nodiscard]] constexpr Vec3 ray_direction(const RayLike auto& ray) {
    return Vec3{
        static_cast<double>(ray.dx),
        static_cast<double>(ray.dy),
        static_cast<double>(ray.dz),
    };
}

[[nodiscard]] constexpr double ray_wavelength(const RayLike auto& ray) {
    return static_cast<double>(ray.wavelength);
}

constexpr void set_ray_origin_mm(RayLike auto& ray, const Vec3 origin_mm) {
    ray.ox = origin_mm.x;
    ray.oy = origin_mm.y;
    ray.oz = origin_mm.z;
}

constexpr void set_ray_direction(RayLike auto& ray, const Vec3 direction) {
    ray.dx = direction.x;
    ray.dy = direction.y;
    ray.dz = direction.z;
}

} // namespace opsys
