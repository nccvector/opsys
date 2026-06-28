#include "osys/osys.hpp"

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

bool traces_plane_parallel_plate() {
    osys::OpticalSystem system{osys::air_medium};
    osys::add_surface(system, osys::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 10.0,
        .medium_after = osys::dense_medium,
        .shape = osys::SagittaSurface{osys::PlaneSagitta{}},
    });
    osys::add_surface(system, osys::OpticalSurface{
        .vertex_z_mm = 5.0,
        .aperture_radius_mm = 10.0,
        .medium_after = osys::air_medium,
        .shape = osys::SagittaSurface{osys::PlaneSagitta{}},
    });

    const osys::Ray ray{
        .origin_mm = {0.0, 0.0, -1.0},
        .direction = {0.09950371902099893, 0.0, 0.9950371902099893},
        .wavelength_nm = 550.0,
    };
    const osys::TraceResult result = osys::trace(system, ray);

    return expect(result.status == osys::TraceStatus::ok, "plane plate should trace")
        && expect(near(result.output_ray.direction.x, ray.direction.x, 1.0e-12), "parallel plate preserves outgoing x angle")
        && expect(near(result.output_ray.direction.z, ray.direction.z, 1.0e-12), "parallel plate preserves outgoing z angle");
}

bool rejects_aperture_miss() {
    osys::OpticalSystem system{osys::air_medium};
    osys::add_surface(system, osys::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 1.0,
        .medium_after = osys::dense_medium,
        .shape = osys::SagittaSurface{osys::PlaneSagitta{}},
    });

    const osys::Ray ray{
        .origin_mm = {2.0, 0.0, -1.0},
        .direction = {0.0, 0.0, 1.0},
        .wavelength_nm = 550.0,
    };
    const osys::TraceResult result = osys::trace(system, ray);

    return expect(result.status == osys::TraceStatus::missed_aperture, "ray outside clear aperture should miss");
}

bool traces_spherical_surface() {
    osys::OpticalSystem system{osys::air_medium};
    osys::add_surface(system, osys::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 10.0,
        .medium_after = osys::dense_medium,
        .shape = osys::SagittaSurface{osys::ConicSagitta{.radius_mm = 50.0}},
    });

    const osys::Ray ray{
        .origin_mm = {0.0, 4.0, -20.0},
        .direction = {0.0, 0.0, 1.0},
        .wavelength_nm = 550.0,
    };
    const osys::TraceResult result = osys::trace(system, ray);

    return expect(result.status == osys::TraceStatus::ok, "spherical surface should trace")
        && expect(result.output_ray.direction.y < 0.0, "positive-radius spherical surface bends marginal ray toward axis");
}

bool traces_paraboloid_surface() {
    osys::OpticalSystem system{osys::air_medium};
    osys::add_surface(system, osys::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 10.0,
        .medium_after = osys::dense_medium,
        .shape = osys::SagittaSurface{osys::ConicSagitta{
            .radius_mm = 50.0,
            .conic_constant = -1.0,
        }},
    });

    const osys::Ray ray{
        .origin_mm = {0.0, 4.0, -20.0},
        .direction = {0.0, 0.0, 1.0},
        .wavelength_nm = 550.0,
    };
    const osys::TraceResult result = osys::trace(system, ray);

    return expect(result.status == osys::TraceStatus::ok, "paraboloid surface should trace")
        && expect(result.output_ray.direction.y < 0.0, "positive-radius paraboloid bends marginal ray toward axis");
}

bool bk7_dispersion_changes_index_with_nm_input() {
    const osys::Medium bk7 = osys::n_bk7_medium;
    const double blue = osys::refractive_index(bk7, 486.1);
    const double red = osys::refractive_index(bk7, 656.3);

    return expect(blue > red, "BK7 should have higher index at shorter wavelength");
}

bool abbe_medium_uses_v_number_dispersion() {
    const osys::Medium medium{osys::AbbeVd{.nd = 1.5168, .vd = 64.1}};
    const double blue = osys::refractive_index(medium, 486.1327);
    const double yellow = osys::refractive_index(medium, 587.5618);
    const double red = osys::refractive_index(medium, 656.2725);

    return expect(blue > yellow, "Abbe medium should bend blue more than yellow")
        && expect(yellow > red, "Abbe medium should bend yellow more than red")
        && expect(near(yellow, 1.5168, 1.0e-6), "Abbe medium should match nd at the d line");
}

bool copies_and_traces_fixed_optical_presets() {
    bool ok = true;
    for (const osys::OpticalPresetId id : osys::optical_preset_ids) {
        const std::span<const osys::OpticalSurfaceSpec> specs = osys::optical_preset_surfaces(id);
        const osys::OpticalSystem system = osys::optical_system(id);
        ok = expect(!specs.empty(), "preset surface table should not be empty") && ok;
        ok = expect(system.surfaces.size() == specs.size(), "preset copy should preserve surface count") && ok;

        const osys::Ray ray{
            .origin_mm = {0.0, 0.0, specs.front().vertex_z_mm - 20.0},
            .direction = {0.0, 0.0, 1.0},
            .wavelength_nm = 550.0,
        };
        const osys::TraceResult result = osys::trace(system, ray);
        ok = expect(result.status == osys::TraceStatus::ok, "center ray should trace through preset") && ok;
    }

    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok = traces_plane_parallel_plate() && ok;
    ok = rejects_aperture_miss() && ok;
    ok = traces_spherical_surface() && ok;
    ok = traces_paraboloid_surface() && ok;
    ok = bk7_dispersion_changes_index_with_nm_input() && ok;
    ok = abbe_medium_uses_v_number_dispersion() && ok;
    ok = copies_and_traces_fixed_optical_presets() && ok;

    return ok ? 0 : 1;
}
