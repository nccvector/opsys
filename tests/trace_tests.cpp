#include "lenses/material.hpp"
#include "lenses/optical_system.hpp"

#include <cmath>
#include <iostream>
#include <string_view>

namespace {

bool expect(bool condition, std::string_view message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        return false;
    }
    return true;
}

bool near(double lhs, double rhs, double tolerance) {
    return std::abs(lhs - rhs) <= tolerance;
}

lenses::MaterialCatalog make_catalog() {
    lenses::MaterialCatalog materials;
    lenses::add_material(materials, "air", lenses::make_air());
    lenses::add_material(materials, "N-BK7", lenses::make_n_bk7());
    lenses::add_material(materials, "dense", lenses::Material{lenses::ConstantIndex{1.5}});
    return materials;
}

bool traces_plane_parallel_plate() {
    lenses::OpticalSystem system{"air"};
    lenses::add_surface(system, lenses::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 10.0,
        .medium_after = "dense",
        .shape = lenses::SagSurface{lenses::PlaneSag{}},
    });
    lenses::add_surface(system, lenses::OpticalSurface{
        .vertex_z_mm = 5.0,
        .aperture_radius_mm = 10.0,
        .medium_after = "air",
        .shape = lenses::SagSurface{lenses::PlaneSag{}},
    });

    const lenses::Ray ray = lenses::make_ray({0.0, 0.0, -1.0}, {0.1, 0.0, 1.0}, 550.0);
    const lenses::TraceResult result = lenses::trace(system, ray, make_catalog());

    return expect(result.status == lenses::TraceStatus::ok, "plane plate should trace")
        && expect(near(result.output_ray.direction.x, ray.direction.x, 1.0e-12), "parallel plate preserves outgoing x angle")
        && expect(near(result.output_ray.direction.z, ray.direction.z, 1.0e-12), "parallel plate preserves outgoing z angle");
}

bool rejects_aperture_miss() {
    lenses::OpticalSystem system{"air"};
    lenses::add_surface(system, lenses::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 1.0,
        .medium_after = "dense",
        .shape = lenses::SagSurface{lenses::PlaneSag{}},
    });

    const lenses::Ray ray = lenses::make_ray({2.0, 0.0, -1.0}, {0.0, 0.0, 1.0}, 550.0);
    const lenses::TraceResult result = lenses::trace(system, ray, make_catalog());

    return expect(result.status == lenses::TraceStatus::missed_aperture, "ray outside clear aperture should miss");
}

bool traces_spherical_surface() {
    lenses::OpticalSystem system{"air"};
    lenses::add_surface(system, lenses::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 10.0,
        .medium_after = "dense",
        .shape = lenses::SagSurface{lenses::SphericalSag{.radius_mm = 50.0}},
    });

    const lenses::Ray ray = lenses::make_ray({0.0, 4.0, -20.0}, {0.0, 0.0, 1.0}, 550.0);
    const lenses::TraceResult result = lenses::trace(system, ray, make_catalog());

    return expect(result.status == lenses::TraceStatus::ok, "spherical surface should trace")
        && expect(result.output_ray.direction.y < 0.0, "positive-radius spherical surface bends marginal ray toward axis");
}

bool bk7_dispersion_changes_index_with_nm_input() {
    const lenses::Material bk7 = lenses::make_n_bk7();
    const double blue = lenses::refractive_index(bk7, 486.1);
    const double red = lenses::refractive_index(bk7, 656.3);

    return expect(blue > red, "BK7 should have higher index at shorter wavelength");
}

} // namespace

int main() {
    bool ok = true;
    ok = traces_plane_parallel_plate() && ok;
    ok = rejects_aperture_miss() && ok;
    ok = traces_spherical_surface() && ok;
    ok = bk7_dispersion_changes_index_with_nm_input() && ok;

    return ok ? 0 : 1;
}
