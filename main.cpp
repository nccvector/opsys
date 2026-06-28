#include "opsys/opsys.hpp"

#include <iostream>

struct DemoRay {
    double ox{};
    double oy{};
    double oz{};
    double dx{};
    double dy{};
    double dz{};
    double wavelength{};
};

int main() {
    opsys::OpticalSystem system = opsys::optical_system(opsys::OpticalPresetId::double_gauss_50mm_f2);

    constexpr DemoRay input{
        .ox = 0.0,
        .oy = 5.0,
        .oz = -20.0,
        .dx = 0.0,
        .dy = 0.0,
        .dz = 1.0,
        .wavelength = 550.0,
    };

    const opsys::TraceResult<DemoRay> result = opsys::trace(system, input);
    std::cout << "trace status: " << static_cast<int>(result.status) << '\n';
    std::cout << "origin mm: "
              << result.output_ray.ox << ", "
              << result.output_ray.oy << ", "
              << result.output_ray.oz << '\n';
    std::cout << "direction: "
              << result.output_ray.dx << ", "
              << result.output_ray.dy << ", "
              << result.output_ray.dz << '\n';

    return result.status == opsys::TraceStatus::ok ? 0 : 1;
}
