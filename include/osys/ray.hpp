#pragma once

#include "osys/geometry.hpp"

namespace osys {

struct Ray {
    Vec3 origin_mm{};
    Vec3 direction{};
    double wavelength_nm{};
};

} // namespace osys
