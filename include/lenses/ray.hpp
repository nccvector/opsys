#pragma once

#include "lenses/geometry.hpp"

namespace lenses {

struct Ray {
    Vec3 origin_mm{};
    Vec3 direction{};
    double wavelength_nm{};
};

[[nodiscard]] inline Ray make_ray(Vec3 origin_mm, Vec3 direction, double wavelength_nm) {
    return Ray{
        .origin_mm = origin_mm,
        .direction = normalize(direction),
        .wavelength_nm = wavelength_nm,
    };
}

} // namespace lenses
