#pragma once

#include "lenses/geometry.hpp"

namespace lenses {

struct Ray {
    Vec3 origin_mm{};
    Vec3 direction{};
    double wavelength_nm{};
};

} // namespace lenses
