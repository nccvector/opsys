#include "lenses/material.hpp"
#include "lenses/optical_system.hpp"

#include <iostream>

int main() {
    lenses::OpticalSystem system{lenses::fixed_medium_catalog.air};
    lenses::add_surface(system, lenses::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 12.5,
        .medium_after = lenses::fixed_medium_catalog.n_bk7,
        .shape = lenses::SagSurface{lenses::SphericalSag{.radius_mm = 50.0}},
    });
    lenses::add_surface(system, lenses::OpticalSurface{
        .vertex_z_mm = 5.0,
        .aperture_radius_mm = 12.5,
        .medium_after = lenses::fixed_medium_catalog.air,
        .shape = lenses::SagSurface{lenses::SphericalSag{.radius_mm = -50.0}},
    });

    const lenses::Ray input{
        .origin_mm = {0.0, 5.0, -20.0},
        .direction = {0.0, 0.0, 1.0},
        .wavelength_nm = 550.0,
    };

    const lenses::TraceResult result = lenses::trace(system, input);
    std::cout << "trace status: " << static_cast<int>(result.status) << '\n';
    std::cout << "origin mm: "
              << result.output_ray.origin_mm.x << ", "
              << result.output_ray.origin_mm.y << ", "
              << result.output_ray.origin_mm.z << '\n';
    std::cout << "direction: "
              << result.output_ray.direction.x << ", "
              << result.output_ray.direction.y << ", "
              << result.output_ray.direction.z << '\n';

    return result.status == lenses::TraceStatus::ok ? 0 : 1;
}
