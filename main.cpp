#include "osys/osys.hpp"

#include <iostream>

int main() {
    osys::OpticalSystem system{osys::air_medium};
    osys::add_surface(system, osys::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 12.5,
        .medium_after = osys::n_bk7_medium,
        .shape = osys::SagittaSurface{osys::ConicSagitta{.radius_mm = 50.0}},
    });
    osys::add_surface(system, osys::OpticalSurface{
        .vertex_z_mm = 5.0,
        .aperture_radius_mm = 12.5,
        .medium_after = osys::air_medium,
        .shape = osys::SagittaSurface{osys::ConicSagitta{.radius_mm = -50.0}},
    });

    constexpr osys::Ray input{
        .origin_mm = {0.0, 5.0, -20.0},
        .direction = {0.0, 0.0, 1.0},
        .wavelength_nm = 550.0,
    };

    const osys::TraceResult result = osys::trace(system, input);
    std::cout << "trace status: " << static_cast<int>(result.status) << '\n';
    std::cout << "origin mm: "
              << result.output_ray.origin_mm.x << ", "
              << result.output_ray.origin_mm.y << ", "
              << result.output_ray.origin_mm.z << '\n';
    std::cout << "direction: "
              << result.output_ray.direction.x << ", "
              << result.output_ray.direction.y << ", "
              << result.output_ray.direction.z << '\n';

    return result.status == osys::TraceStatus::ok ? 0 : 1;
}
