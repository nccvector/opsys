# osys

`osys` is a small header-only C++23 optical-system library. It traces rays through fixed optical systems made from value-type surfaces, mediums, and sagitta models. The current implementation favors compile-time-known data, pure free functions, copied medium structs, and analytic plane/conic intersections instead of runtime lookup tables or inheritance-heavy object models.

The repository also includes a raylib editor for inspecting and tweaking optical prescriptions. The editor can export deterministic PNG visualizations for every fixed preset.

## What is included

- Header-only `osys` library under `include/osys`.
- Value-based `Medium`, `OpticalSurface`, `OpticalSystem`, and `Ray` types.
- `PlaneSagitta` and `ConicSagitta` surface models.
- Fixed compile-time optical preset tables in `include/osys/presets.hpp`.
- Spectral ray tracing with wavelength-dependent medium indices.
- Optional raylib editor target, `osys_editor`.
- Tests for tracing, preset copying, and medium dispersion.

## Build

```sh
cmake -S . -B cmake-build-debug
cmake --build cmake-build-debug
ctest --test-dir cmake-build-debug --output-on-failure
```

If raylib is available, CMake also builds the editor:

```sh
./cmake-build-debug/osys_editor
```

Export the preset gallery:

```sh
./cmake-build-debug/osys_editor --export-presets docs/presets
```

## Minimal use

```cpp
#include <osys/osys.hpp>

int main() {
    const osys::OpticalSystem system = osys::optical_system(osys::OpticalPresetId::double_gauss_50mm_f2);
    const osys::Ray ray{
        .origin_mm = {0.0, 2.0, -40.0},
        .direction = {0.0, 0.0, 1.0},
        .wavelength_nm = 550.0,
    };

    const osys::TraceResult result = osys::trace(system, ray);
    return result.status == osys::TraceStatus::ok ? 0 : 1;
}
```

## Fixed preset gallery

| Preset | Visualization |
| --- | --- |
| Double-Gauss 50mm f/2 | <img src="docs/presets/double-gauss-50mm-f-2.png" width="420" alt="Double-Gauss 50mm f/2 preset visualization"> |
| Nakamura wide-angle 22mm | <img src="docs/presets/nakamura-wide-angle-22mm.png" width="420" alt="Nakamura wide-angle 22mm preset visualization"> |
| Muller fisheye 10mm f/4 | <img src="docs/presets/muller-fisheye-10mm-f-4.png" width="420" alt="Muller fisheye 10mm f/4 preset visualization"> |
| Sigler telephoto 250mm f/5.6 | <img src="docs/presets/sigler-telephoto-250mm-f-5-6.png" width="420" alt="Sigler telephoto 250mm f/5.6 preset visualization"> |
| Nikon Nikkor-S Auto 50mm f/1.4 | <img src="docs/presets/nikon-nikkor-s-auto-50mm-f-1-4.png" width="420" alt="Nikon Nikkor-S Auto 50mm f/1.4 preset visualization"> |
| Nikon AI Nikkor 28mm f/2.8S | <img src="docs/presets/nikon-ai-nikkor-28mm-f-2-8s.png" width="420" alt="Nikon AI Nikkor 28mm f/2.8S preset visualization"> |
| Nikon AI AF Fisheye-Nikkor 16mm f/2.8D | <img src="docs/presets/nikon-ai-af-fisheye-nikkor-16mm-f-2-8d.png" width="420" alt="Nikon AI AF Fisheye-Nikkor 16mm f/2.8D preset visualization"> |
| Nikon AI Nikkor 180mm f/2.8 ED | <img src="docs/presets/nikon-ai-nikkor-180mm-f-2-8-ed.png" width="420" alt="Nikon AI Nikkor 180mm f/2.8 ED preset visualization"> |

## Design constraints

`osys` deliberately avoids runtime OOP. Mediums and surfaces are small copied values, preset prescriptions are fixed tables, and tracing reads the optical system directly. Dynamic catalogs and runtime prescription loading can be added later as a separate layer without changing the core value-oriented tracing path.
