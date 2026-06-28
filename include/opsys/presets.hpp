#pragma once

#include "opsys/optical_system.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <string_view>

namespace opsys {

struct OpticalSurfaceSpec {
    double vertex_z_mm{};
    double aperture_radius_mm{};
    Medium medium_after{};
    SagittaSurface shape{};
    std::string_view medium_name_after{};
};

enum class OpticalPresetId {
    double_gauss_50mm_f2,
    wide_angle_22mm,
    fisheye_10mm,
    telephoto_250mm,
    nikkor_50mm_f1_4,
    nikkor_28mm_f2_8s,
    nikkor_16mm_f2_8d,
    nikkor_180mm_f2_8_ed,
};

struct OpticalPresetInfo {
    OpticalPresetId id{};
    std::string_view name{};
    std::string_view source{};
    double focal_length_mm{};
    double f_number{};
};

[[nodiscard]] inline constexpr Medium nd_medium(const double n) {
    return Medium{ConstantIndex{n}};
}

[[nodiscard]] inline constexpr Medium abbe_medium(const double nd, const double vd) {
    return Medium{AbbeVd{.nd = nd, .vd = vd}};
}

[[nodiscard]] inline constexpr SagittaSurface spherical_sagitta(const double radius_mm) {
    return SagittaSurface{ConicSagitta{.radius_mm = radius_mm, .conic_constant = 0.0}};
}

[[nodiscard]] inline constexpr OpticalSurfaceSpec spherical_surface(
    const double vertex_z_mm,
    const double aperture_radius_mm,
    const Medium medium_after,
    const std::string_view medium_name_after,
    const double radius_mm) {
    return OpticalSurfaceSpec{
        .vertex_z_mm = vertex_z_mm,
        .aperture_radius_mm = aperture_radius_mm,
        .medium_after = medium_after,
        .shape = spherical_sagitta(radius_mm),
        .medium_name_after = medium_name_after,
    };
}

[[nodiscard]] inline constexpr OpticalSurfaceSpec plane_surface(
    const double vertex_z_mm,
    const double aperture_radius_mm,
    const Medium medium_after,
    const std::string_view medium_name_after) {
    return OpticalSurfaceSpec{
        .vertex_z_mm = vertex_z_mm,
        .aperture_radius_mm = aperture_radius_mm,
        .medium_after = medium_after,
        .shape = SagittaSurface{PlaneSagitta{}},
        .medium_name_after = medium_name_after,
    };
}

[[nodiscard]] inline constexpr OpticalSurface optical_surface(const OpticalSurfaceSpec& spec) {
    return OpticalSurface{
        .vertex_z_mm = spec.vertex_z_mm,
        .aperture_radius_mm = spec.aperture_radius_mm,
        .medium_after = spec.medium_after,
        .shape = spec.shape,
    };
}

inline constexpr std::array<OpticalSurfaceSpec, 11> double_gauss_50mm_f2_surfaces{
    spherical_surface(0.0, 12.6, nd_medium(1.67), "n 1.670", 29.475),
    spherical_surface(3.76, 12.6, air_medium, "Air", 84.83),
    spherical_surface(3.88, 11.5, nd_medium(1.67), "n 1.670", 19.275),
    spherical_surface(7.905, 11.5, nd_medium(1.699), "n 1.699", 40.77),
    spherical_surface(11.18, 9.0, air_medium, "Air", 12.75),
    plane_surface(16.885, 8.55, air_medium, "Air"),
    spherical_surface(21.385, 8.5, nd_medium(1.603), "n 1.603", -14.495),
    spherical_surface(22.565, 10.0, nd_medium(1.658), "n 1.658", 40.77),
    spherical_surface(28.63, 10.0, air_medium, "Air", -20.385),
    spherical_surface(28.82, 10.0, nd_medium(1.717), "n 1.717", 437.065),
    spherical_surface(32.04, 10.0, air_medium, "Air", -39.73),
};

inline constexpr std::array<OpticalSurfaceSpec, 13> wide_angle_22mm_surfaces{
    spherical_surface(0.0, 11.858, nd_medium(1.54), "n 1.540", 35.98738),
    spherical_surface(1.21638, 8.998, air_medium, "Air", 11.69718),
    spherical_surface(11.21208, 6.182, nd_medium(1.772), "n 1.772", 13.08714),
    spherical_surface(16.3383, 4.906, nd_medium(1.617), "n 1.617", -22.63294),
    spherical_surface(18.10754, 4.576, air_medium, "Air", 71.05802),
    plane_surface(18.92594, 4.378, air_medium, "Air"),
    spherical_surface(21.2036, 4.092, nd_medium(1.617), "n 1.617", -9.58584),
    spherical_surface(23.63614, 4.576, air_medium, "Air", -11.28864),
    spherical_surface(23.7512, 5.324, nd_medium(1.713), "n 1.713", -166.7765),
    spherical_surface(26.84726, 5.72, nd_medium(1.805), "n 1.805", -7.5911),
    spherical_surface(28.17408, 6.138, air_medium, "Air", -16.7662),
    spherical_surface(32.15476, 6.71, nd_medium(1.617), "n 1.617", -7.70286),
    spherical_surface(33.37114, 8.998, air_medium, "Air", -11.97328),
};

inline constexpr std::array<OpticalSurfaceSpec, 12> fisheye_10mm_surfaces{
    spherical_surface(0.0, 15.17, nd_medium(1.62), "n 1.620", 30.2249),
    spherical_surface(0.8335, 10.34, air_medium, "Air", 11.3931),
    spherical_surface(8.2471, 8.9, nd_medium(1.639), "n 1.639", 75.2019),
    spherical_surface(9.3125, 6.71, air_medium, "Air", 8.3349),
    spherical_surface(20.4674, 4.51, nd_medium(1.654), "n 1.654", 9.5882),
    spherical_surface(22.4728, 4.07, air_medium, "Air", 43.8677),
    plane_surface(27.8623, 3.04, air_medium, "Air"),
    spherical_surface(29.2786, 2.98, nd_medium(1.517), "n 1.517", 29.4541),
    spherical_surface(31.472, 2.92, nd_medium(1.805), "n 1.805", -5.2265),
    spherical_surface(32.4434, 2.98, air_medium, "Air", -14.2884),
    spherical_surface(32.5061, 2.98, nd_medium(1.673), "n 1.673", -22.3726),
    spherical_surface(33.4461, 3.26, air_medium, "Air", -15.0404),
};

inline constexpr std::array<OpticalSurfaceSpec, 7> telephoto_250mm_surfaces{
    spherical_surface(0.0, 23.75, nd_medium(1.529), "n 1.529", 54.6275),
    spherical_surface(12.52, 22.25, nd_medium(1.599), "n 1.599", -86.365),
    spherical_surface(16.275, 20.75, air_medium, "Air", 271.7625),
    plane_surface(19.0925, 20.25, air_medium, "Air"),
    spherical_surface(86.505, 15.75, nd_medium(1.613), "n 1.613", -32.13),
    spherical_surface(90.26, 16.75, nd_medium(1.603), "n 1.603", 49.5325),
    spherical_surface(102.78, 18.5, air_medium, "Air", -50.945),
};

inline constexpr std::array<OpticalSurfaceSpec, 13> nikkor_50mm_f1_4_surfaces{
    spherical_surface(0.0, 19.24, abbe_medium(1.66446, 35.9), "nD1.664/V35.9", 38.275),
    spherical_surface(6.3, 19.24, air_medium, "Air", 125.265),
    spherical_surface(6.395, 17.26, abbe_medium(1.69350, 53.4), "nD1.694/V53.4", 28.68),
    spherical_surface(14.825, 17.26, abbe_medium(1.64831, 33.8), "nD1.648/V33.8", 620.205),
    spherical_surface(16.47, 11.97, air_medium, "Air", 16.2),
    plane_surface(24.29, 11.6895, air_medium, "Air"),
    spherical_surface(32.75, 11.475, abbe_medium(1.69895, 30.1), "nD1.699/V30.1", -14.46),
    spherical_surface(34.785, 15.36, abbe_medium(1.71300, 53.9), "nD1.713/V53.9", -193.8),
    spherical_surface(43.02, 15.36, air_medium, "Air", -22.675),
    spherical_surface(43.215, 16.375, abbe_medium(1.67790, 55.5), "nD1.678/V55.5", -116.28),
    spherical_surface(47.38, 16.375, air_medium, "Air", -31.785),
    spherical_surface(47.475, 17.26, abbe_medium(1.71300, 53.9), "nD1.713/V53.9", 77.52),
    spherical_surface(50.575, 17.26, air_medium, "Air", -843.87),
};

inline constexpr std::array<OpticalSurfaceSpec, 17> nikkor_28mm_f2_8s_surfaces{
    spherical_surface(0.0, 20.535, abbe_medium(1.67003, 47.07), "nD1.670/V47.1", 60.5),
    spherical_surface(3.5, 20.535, air_medium, "Air", 118.3),
    spherical_surface(3.6, 17.335, abbe_medium(1.51680, 64.10), "nD1.517/V64.1", 36.9),
    spherical_surface(5.1, 13.0, air_medium, "Air", 15.4),
    spherical_surface(9.1, 14.005, abbe_medium(1.51680, 64.10), "nD1.517/V64.1", 28.2),
    spherical_surface(10.6, 12.02, air_medium, "Air", 15.13),
    spherical_surface(25.0, 7.75, abbe_medium(1.66755, 41.96), "nD1.668/V42.0", 28.6),
    plane_surface(29.4, 7.75, air_medium, "Air"),
    plane_surface(30.75, 10.275, abbe_medium(1.62041, 60.14), "nD1.620/V60.1"),
    spherical_surface(34.75, 10.275, air_medium, "Air", -33.97),
    plane_surface(35.75, 6.915, air_medium, "Air"),
    spherical_surface(40.05, 7.94, abbe_medium(1.75520, 27.61), "nD1.755/V27.6", -20.8),
    spherical_surface(42.6, 7.0, air_medium, "Air", 53.18),
    spherical_surface(43.6, 7.0, abbe_medium(1.67025, 57.53), "nD1.670/V57.5", -46.001),
    spherical_surface(46.7, 9.07, air_medium, "Air", -17.9),
    spherical_surface(46.8, 10.555, abbe_medium(1.67025, 57.53), "nD1.670/V57.5", 109.097),
    spherical_surface(50.3, 10.555, air_medium, "Air", -34.8),
};

inline constexpr std::array<OpticalSurfaceSpec, 16> nikkor_16mm_f2_8d_surfaces{
    spherical_surface(0.0, 26.305, abbe_medium(1.64000, 60.0), "nD1.640/V60.0", 69.258),
    spherical_surface(1.7, 14.65, air_medium, "Air", 15.165),
    spherical_surface(11.95, 15.135, abbe_medium(1.62041, 60.1), "nD1.620/V60.1", 130.305),
    spherical_surface(13.65, 11.41, air_medium, "Air", 14.794),
    spherical_surface(17.65, 10.135, abbe_medium(1.62588, 35.7), "nD1.626/V35.7", 32.412),
    spherical_surface(25.65, 10.135, abbe_medium(1.79668, 45.4), "nD1.797/V45.4", -15.298),
    spherical_surface(34.65, 7.24, air_medium, "Air", -38.75),
    plane_surface(41.7, 7.1785, air_medium, "Air"),
    spherical_surface(45.63, 7.7, abbe_medium(1.58144, 40.8), "nD1.581/V40.8", -81.004),
    spherical_surface(46.83, 8.305, abbe_medium(1.51823, 58.9), "nD1.518/V58.9", 24.002),
    spherical_surface(51.83, 8.305, air_medium, "Air", -20.608),
    spherical_surface(52.72, 9.815, abbe_medium(1.51860, 70.0), "nD1.519/V70.0", 46.733),
    spherical_surface(58.22, 9.815, abbe_medium(1.78470, 26.1), "nD1.785/V26.1", -19.219),
    spherical_surface(59.52, 10.83, air_medium, "Air", -40.706),
    plane_surface(60.37, 11.605, abbe_medium(1.51680, 64.1), "nD1.517/V64.1"),
    plane_surface(61.57, 11.605, air_medium, "Air"),
};

inline constexpr std::array<OpticalSurfaceSpec, 11> nikkor_180mm_f2_8_ed_surfaces{
    spherical_surface(0.0, 32.47, abbe_medium(1.50032, 81.9), "nD1.500/V81.9", 99.022),
    spherical_surface(11.5, 32.47, air_medium, "Air", -140.839),
    spherical_surface(13.6, 32.47, abbe_medium(1.74950, 35.0), "nD1.750/V35.0", -138.056),
    spherical_surface(17.3, 32.47, air_medium, "Air", 373.0),
    spherical_surface(23.6, 30.5, abbe_medium(1.65844, 50.8), "nD1.658/V50.8", 77.774),
    spherical_surface(32.8, 30.5, air_medium, "Air", 239.999),
    plane_surface(90.9, 15.7515, air_medium, "Air"),
    spherical_surface(107.3, 15.68, abbe_medium(1.51454, 54.6), "nD1.515/V54.6", -35.5),
    spherical_surface(109.1, 15.68, air_medium, "Air", -550.001),
    spherical_surface(109.6, 15.68, abbe_medium(1.79668, 45.4), "nD1.797/V45.4", 220.0),
    spherical_surface(114.6, 15.68, air_medium, "Air", -162.193),
};

inline constexpr std::array<OpticalPresetId, 8> optical_preset_ids{
    OpticalPresetId::double_gauss_50mm_f2,
    OpticalPresetId::wide_angle_22mm,
    OpticalPresetId::fisheye_10mm,
    OpticalPresetId::telephoto_250mm,
    OpticalPresetId::nikkor_50mm_f1_4,
    OpticalPresetId::nikkor_28mm_f2_8s,
    OpticalPresetId::nikkor_16mm_f2_8d,
    OpticalPresetId::nikkor_180mm_f2_8_ed,
};

inline constexpr std::array<OpticalPresetInfo, optical_preset_ids.size()> optical_preset_infos{
    OpticalPresetInfo{OpticalPresetId::double_gauss_50mm_f2, "Double-Gauss 50mm f/2", "PBRT v4 scenes / US 2,673,491 / Modern Lens Design", 50.0, 2.0},
    OpticalPresetInfo{OpticalPresetId::wide_angle_22mm, "Nakamura wide-angle 22mm", "PBRT v4 scenes / Modern Lens Design p.360", 22.0, 2.5},
    OpticalPresetInfo{OpticalPresetId::fisheye_10mm, "Muller fisheye 10mm f/4", "PBRT v4 scenes / Modern Lens Design p.164", 10.0, 4.0},
    OpticalPresetInfo{OpticalPresetId::telephoto_250mm, "Sigler telephoto 250mm f/5.6", "PBRT v4 scenes / Modern Lens Design p.175", 250.0, 5.6},
    OpticalPresetInfo{OpticalPresetId::nikkor_50mm_f1_4, "Nikon Nikkor-S Auto 50mm f/1.4", "Photons to Photos Optical Bench / US 3,560,079 Example 1", 50.0, 1.4},
    OpticalPresetInfo{OpticalPresetId::nikkor_28mm_f2_8s, "Nikon AI Nikkor 28mm f/2.8S", "Photons to Photos Optical Bench / US 5,917,663 Example 2", 28.61, 2.89},
    OpticalPresetInfo{OpticalPresetId::nikkor_16mm_f2_8d, "Nikon AI AF Fisheye-Nikkor 16mm f/2.8D", "Photons to Photos Optical Bench / US 5,434,713 Example 2", 15.6686, 2.8},
    OpticalPresetInfo{OpticalPresetId::nikkor_180mm_f2_8_ed, "Nikon AI Nikkor 180mm f/2.8 ED", "Photons to Photos Optical Bench / US 4,338,001 Example 2", 180.0, 2.8},
};

[[nodiscard]] inline constexpr std::size_t optical_preset_index(const OpticalPresetId id) {
    for (std::size_t i = 0; i < optical_preset_ids.size(); ++i) {
        if (optical_preset_ids[i] == id) {
            return i;
        }
    }

    return 0;
}

[[nodiscard]] inline constexpr OpticalPresetId next_optical_preset_id(const OpticalPresetId id) {
    const std::size_t index = optical_preset_index(id);
    return optical_preset_ids[(index + 1) % optical_preset_ids.size()];
}

[[nodiscard]] inline constexpr OpticalPresetId previous_optical_preset_id(const OpticalPresetId id) {
    const std::size_t index = optical_preset_index(id);
    return optical_preset_ids[(index + optical_preset_ids.size() - 1) % optical_preset_ids.size()];
}

[[nodiscard]] inline constexpr const OpticalPresetInfo& optical_preset_info(const OpticalPresetId id) {
    return optical_preset_infos[optical_preset_index(id)];
}

[[nodiscard]] inline constexpr std::span<const OpticalSurfaceSpec> optical_preset_surfaces(const OpticalPresetId id) {
    switch (id) {
        case OpticalPresetId::double_gauss_50mm_f2:
            return std::span<const OpticalSurfaceSpec>{double_gauss_50mm_f2_surfaces};
        case OpticalPresetId::wide_angle_22mm:
            return std::span<const OpticalSurfaceSpec>{wide_angle_22mm_surfaces};
        case OpticalPresetId::fisheye_10mm:
            return std::span<const OpticalSurfaceSpec>{fisheye_10mm_surfaces};
        case OpticalPresetId::telephoto_250mm:
            return std::span<const OpticalSurfaceSpec>{telephoto_250mm_surfaces};
        case OpticalPresetId::nikkor_50mm_f1_4:
            return std::span<const OpticalSurfaceSpec>{nikkor_50mm_f1_4_surfaces};
        case OpticalPresetId::nikkor_28mm_f2_8s:
            return std::span<const OpticalSurfaceSpec>{nikkor_28mm_f2_8s_surfaces};
        case OpticalPresetId::nikkor_16mm_f2_8d:
            return std::span<const OpticalSurfaceSpec>{nikkor_16mm_f2_8d_surfaces};
        case OpticalPresetId::nikkor_180mm_f2_8_ed:
            return std::span<const OpticalSurfaceSpec>{nikkor_180mm_f2_8_ed_surfaces};
    }

    return std::span<const OpticalSurfaceSpec>{double_gauss_50mm_f2_surfaces};
}

[[nodiscard]] inline OpticalSystem optical_system(const std::span<const OpticalSurfaceSpec> surfaces) {
    OpticalSystem system{air_medium};
    system.surfaces.reserve(surfaces.size());
    for (const OpticalSurfaceSpec& surface : surfaces) {
        add_surface(system, optical_surface(surface));
    }
    return system;
}

[[nodiscard]] inline OpticalSystem optical_system(const OpticalPresetId id) {
    return optical_system(optical_preset_surfaces(id));
}

} // namespace opsys
