#include "lenses/material.hpp"
#include "lenses/optical_system.hpp"

#include <iostream>

int main() {
    lenses::MaterialCatalog materials;
    lenses::add_material(materials, "air", lenses::make_air());
    lenses::add_material(materials, "N-BK7", lenses::make_n_bk7());

    lenses::OpticalSystem system{"air"};
    lenses::add_surface(system, lenses::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 12.5,
        .medium_after = "N-BK7",
        .shape = lenses::SagSurface{lenses::SphericalSag{.radius_mm = 50.0}},
    });
    lenses::add_surface(system, lenses::OpticalSurface{
        .vertex_z_mm = 5.0,
        .aperture_radius_mm = 12.5,
        .medium_after = "air",
        .shape = lenses::SagSurface{lenses::SphericalSag{.radius_mm = -50.0}},
    });

    const lenses::Ray input = lenses::make_ray(
        lenses::Vec3{0.0, 5.0, -20.0},
        lenses::Vec3{0.0, 0.0, 1.0},
        550.0);

    const lenses::TraceResult result = lenses::trace(system, input, materials);
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
