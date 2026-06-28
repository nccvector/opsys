#include "opsys/opsys.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <optional>
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

bool near_point(const std::array<double, 3> lhs, const std::array<double, 3> rhs, const double tolerance) {
    return near(lhs[0], rhs[0], tolerance)
        && near(lhs[1], rhs[1], tolerance)
        && near(lhs[2], rhs[2], tolerance);
}

[[nodiscard]] std::array<double, 3> normalized_direction(const opsys::RayLike auto& ray) {
    const double length = std::sqrt(ray.dx * ray.dx + ray.dy * ray.dy + ray.dz * ray.dz);
    return {ray.dx / length, ray.dy / length, ray.dz / length};
}

struct TestRay {
    double ox{};
    double oy{};
    double oz{};
    double dx{};
    double dy{};
    double dz{};
    double wavelength{};
};

std::optional<std::array<double, 3>> point_at_z(const opsys::RayLike auto& ray, const double z_mm) {
    if (near(ray.dz, 0.0, 1.0e-14)) {
        return std::nullopt;
    }

    const double t = (z_mm - ray.oz) / ray.dz;
    if (!std::isfinite(t)) {
        return std::nullopt;
    }

    return std::array<double, 3>{
        ray.ox + t * ray.dx,
        ray.oy + t * ray.dy,
        ray.oz + t * ray.dz,
    };
}

bool traces_plane_parallel_plate() {
    opsys::OpticalSystem system{opsys::air_medium};
    opsys::add_surface(system, opsys::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 10.0,
        .medium_after = opsys::dense_medium,
        .shape = opsys::SagittaSurface{opsys::PlaneSagitta{}},
    });
    opsys::add_surface(system, opsys::OpticalSurface{
        .vertex_z_mm = 5.0,
        .aperture_radius_mm = 10.0,
        .medium_after = opsys::air_medium,
        .shape = opsys::SagittaSurface{opsys::PlaneSagitta{}},
    });

    const TestRay ray{
        .ox = 0.0,
        .oy = 0.0,
        .oz = -1.0,
        .dx = 0.09950371902099893,
        .dy = 0.0,
        .dz = 0.9950371902099893,
        .wavelength = 550.0,
    };
    const opsys::TraceResult<TestRay> result = opsys::trace(system, ray);

    return expect(result.status == opsys::TraceStatus::ok, "plane plate should trace")
        && expect(near(result.output_ray.dx, ray.dx, 1.0e-12), "parallel plate preserves outgoing x angle")
        && expect(near(result.output_ray.dz, ray.dz, 1.0e-12), "parallel plate preserves outgoing z angle");
}

bool rejects_aperture_miss() {
    opsys::OpticalSystem system{opsys::air_medium};
    opsys::add_surface(system, opsys::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 1.0,
        .medium_after = opsys::dense_medium,
        .shape = opsys::SagittaSurface{opsys::PlaneSagitta{}},
    });

    const TestRay ray{
        .ox = 2.0,
        .oy = 0.0,
        .oz = -1.0,
        .dx = 0.0,
        .dy = 0.0,
        .dz = 1.0,
        .wavelength = 550.0,
    };
    const opsys::TraceResult<TestRay> result = opsys::trace(system, ray);

    return expect(result.status == opsys::TraceStatus::missed_aperture, "ray outside clear aperture should miss");
}

bool traces_spherical_surface() {
    opsys::OpticalSystem system{opsys::air_medium};
    opsys::add_surface(system, opsys::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 10.0,
        .medium_after = opsys::dense_medium,
        .shape = opsys::SagittaSurface{opsys::ConicSagitta{.radius_mm = 50.0}},
    });

    const TestRay ray{
        .ox = 0.0,
        .oy = 4.0,
        .oz = -20.0,
        .dx = 0.0,
        .dy = 0.0,
        .dz = 1.0,
        .wavelength = 550.0,
    };
    const opsys::TraceResult<TestRay> result = opsys::trace(system, ray);

    return expect(result.status == opsys::TraceStatus::ok, "spherical surface should trace")
        && expect(result.output_ray.dy < 0.0, "positive-radius spherical surface bends marginal ray toward axis");
}

bool traces_paraboloid_surface() {
    opsys::OpticalSystem system{opsys::air_medium};
    opsys::add_surface(system, opsys::OpticalSurface{
        .vertex_z_mm = 0.0,
        .aperture_radius_mm = 10.0,
        .medium_after = opsys::dense_medium,
        .shape = opsys::SagittaSurface{opsys::ConicSagitta{
            .radius_mm = 50.0,
            .conic_constant = -1.0,
        }},
    });

    const TestRay ray{
        .ox = 0.0,
        .oy = 4.0,
        .oz = -20.0,
        .dx = 0.0,
        .dy = 0.0,
        .dz = 1.0,
        .wavelength = 550.0,
    };
    const opsys::TraceResult<TestRay> result = opsys::trace(system, ray);

    return expect(result.status == opsys::TraceStatus::ok, "paraboloid surface should trace")
        && expect(result.output_ray.dy < 0.0, "positive-radius paraboloid bends marginal ray toward axis");
}

bool bk7_dispersion_changes_index_with_nm_input() {
    const opsys::Medium bk7 = opsys::n_bk7_medium;
    const double blue = opsys::refractive_index(bk7, 486.1);
    const double red = opsys::refractive_index(bk7, 656.3);

    return expect(blue > red, "BK7 should have higher index at shorter wavelength");
}

bool abbe_medium_uses_v_number_dispersion() {
    const opsys::Medium medium{opsys::AbbeVd{.nd = 1.5168, .vd = 64.1}};
    const double blue = opsys::refractive_index(medium, 486.1327);
    const double yellow = opsys::refractive_index(medium, 587.5618);
    const double red = opsys::refractive_index(medium, 656.2725);

    return expect(blue > yellow, "Abbe medium should bend blue more than yellow")
        && expect(yellow > red, "Abbe medium should bend yellow more than red")
        && expect(near(yellow, 1.5168, 1.0e-6), "Abbe medium should match nd at the d line");
}

bool copies_and_traces_fixed_optical_presets() {
    bool ok = true;
    for (const opsys::OpticalPresetId id : opsys::optical_preset_ids) {
        const std::span<const opsys::OpticalSurfaceSpec> specs = opsys::optical_preset_surfaces(id);
        const opsys::OpticalSystem system = opsys::optical_system(id);
        ok = expect(!specs.empty(), "preset surface table should not be empty") && ok;
        ok = expect(system.surfaces.size() == specs.size(), "preset copy should preserve surface count") && ok;

        const TestRay ray{
            .ox = 0.0,
            .oy = 0.0,
            .oz = specs.front().vertex_z_mm - 20.0,
            .dx = 0.0,
            .dy = 0.0,
            .dz = 1.0,
            .wavelength = 550.0,
        };
        const opsys::TraceResult<TestRay> result = opsys::trace(system, ray);
        ok = expect(result.status == opsys::TraceStatus::ok, "center ray should trace through preset") && ok;
    }

    return ok;
}

bool reverse_trace_round_trips_fixed_optical_presets() {
    bool ok = true;
    constexpr std::array<double, 3> wavelengths_nm{486.1327, 587.5618, 656.2725};

    for (const opsys::OpticalPresetId id : opsys::optical_preset_ids) {
        const std::span<const opsys::OpticalSurfaceSpec> specs = opsys::optical_preset_surfaces(id);
        const opsys::OpticalSystem system = opsys::optical_system(id);
        const double world_z_mm = specs.front().vertex_z_mm - 20.0;
        const double image_z_mm = specs.back().vertex_z_mm + 20.0;
        const double off_axis_y_mm = 0.04 * specs.front().aperture_radius_mm;
        const std::array<double, 2> ray_heights_mm{0.0, off_axis_y_mm};

        for (const double wavelength_nm : wavelengths_nm) {
            for (const double ray_height_mm : ray_heights_mm) {
                const TestRay forward_input{
                    .ox = 0.0,
                    .oy = ray_height_mm,
                    .oz = world_z_mm,
                    .dx = 0.0,
                    .dy = 0.0,
                    .dz = 1.0,
                    .wavelength = wavelength_nm,
                };
                const opsys::TraceResult<TestRay> forward = opsys::trace(system, forward_input);
                ok = expect(forward.status == opsys::TraceStatus::ok, "forward preset ray should trace") && ok;

                const std::optional<std::array<double, 3>> image_point = point_at_z(forward.output_ray, image_z_mm);
                ok = expect(image_point.has_value(), "forward output ray should reach image-side plane") && ok;
                if (forward.status != opsys::TraceStatus::ok || !image_point.has_value()) {
                    continue;
                }

                const TestRay reverse_input{
                    .ox = (*image_point)[0],
                    .oy = (*image_point)[1],
                    .oz = (*image_point)[2],
                    .dx = -forward.output_ray.dx,
                    .dy = -forward.output_ray.dy,
                    .dz = -forward.output_ray.dz,
                    .wavelength = wavelength_nm,
                };
                const opsys::TraceResult<TestRay> reverse = opsys::reverse_trace(system, reverse_input);
                ok = expect(reverse.status == opsys::TraceStatus::ok, "reverse preset ray should trace") && ok;

                const std::optional<std::array<double, 3>> returned_world_point = point_at_z(reverse.output_ray, world_z_mm);
                ok = expect(returned_world_point.has_value(), "reverse output ray should reach world-side plane") && ok;
                if (reverse.status != opsys::TraceStatus::ok || !returned_world_point.has_value()) {
                    continue;
                }

                ok = expect(
                    near_point(*returned_world_point, {forward_input.ox, forward_input.oy, forward_input.oz}, 1.0e-6),
                    "reverse ray should return to the forward input point") && ok;
                ok = expect(
                    near_point(normalized_direction(reverse.output_ray), {
                        -normalized_direction(forward_input)[0],
                        -normalized_direction(forward_input)[1],
                        -normalized_direction(forward_input)[2],
                    }, 1.0e-10),
                    "reverse ray should leave opposite the forward input direction") && ok;
            }
        }
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
    ok = reverse_trace_round_trips_fixed_optical_presets() && ok;

    return ok ? 0 : 1;
}
