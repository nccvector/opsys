#pragma once

#include <concepts>

namespace opsys {

template <typename Ray>
concept SpectralRay = requires(Ray ray, const Ray const_ray, double value) {
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

[[nodiscard]] constexpr double ray_wavelength(const SpectralRay auto& ray) {
    return static_cast<double>(ray.wavelength);
}

} // namespace opsys
