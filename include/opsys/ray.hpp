#pragma once

#include <concepts>
#include <type_traits>
#include <utility>

namespace opsys {

namespace detail {

template <typename Ray>
using ray_scalar_t = std::common_type_t<
    std::remove_cvref_t<decltype(std::declval<const Ray&>().ox)>,
    std::remove_cvref_t<decltype(std::declval<const Ray&>().oy)>,
    std::remove_cvref_t<decltype(std::declval<const Ray&>().oz)>,
    std::remove_cvref_t<decltype(std::declval<const Ray&>().dx)>,
    std::remove_cvref_t<decltype(std::declval<const Ray&>().dy)>,
    std::remove_cvref_t<decltype(std::declval<const Ray&>().dz)>,
    std::remove_cvref_t<decltype(std::declval<const Ray&>().wavelength)>>;

} // namespace detail

template <typename Ray>
using ray_scalar_t = detail::ray_scalar_t<Ray>;

template <typename Ray>
concept SpectralRay = requires(Ray ray, const Ray const_ray, ray_scalar_t<Ray> value) {
    requires std::floating_point<ray_scalar_t<Ray>>;
    { const_ray.ox } -> std::convertible_to<ray_scalar_t<Ray>>;
    { const_ray.oy } -> std::convertible_to<ray_scalar_t<Ray>>;
    { const_ray.oz } -> std::convertible_to<ray_scalar_t<Ray>>;
    { const_ray.dx } -> std::convertible_to<ray_scalar_t<Ray>>;
    { const_ray.dy } -> std::convertible_to<ray_scalar_t<Ray>>;
    { const_ray.dz } -> std::convertible_to<ray_scalar_t<Ray>>;
    { const_ray.wavelength } -> std::convertible_to<ray_scalar_t<Ray>>;

    ray.ox = value;
    ray.oy = value;
    ray.oz = value;
    ray.dx = value;
    ray.dy = value;
    ray.dz = value;
    ray.wavelength = value;
};

template <SpectralRay Ray>
[[nodiscard]] constexpr ray_scalar_t<Ray> ray_wavelength(const Ray& ray) {
    return static_cast<ray_scalar_t<Ray>>(ray.wavelength);
}

} // namespace opsys
