#include "osys/osys.hpp"

#include <raylib.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace {

constexpr double k_pi = 3.14159265358979323846;
constexpr double k_ray_epsilon_mm = 1.0e-7;
constexpr int k_window_width = 1280;
constexpr int k_window_height = 800;
constexpr int k_min_window_width = 960;
constexpr int k_min_window_height = 640;

enum class SurfaceKind {
    plane,
    conic,
};

enum class RayPreset {
    single_offset,
    parallel_bundle,
    angled_center,
    edge_pair,
    point_source,
};

struct SurfaceEdit {
    bool expanded{true};
    double vertex_z_mm{};
    double aperture_radius_mm{12.5};
    osys::Medium medium_after{osys::n_bk7_medium};
    std::string medium_name_after{"N-BK7"};
    int medium_catalog_index{1};
    SurfaceKind kind{SurfaceKind::conic};
    double radius_mm{50.0};
    double conic_constant{};
};

struct RayControls {
    RayPreset preset{RayPreset::parallel_bundle};
    bool spectral{true};
    double wavelength_nm{550.0};
    int spectral_samples{13};
    int ray_count{7};
    int point_samples{9};
    double single_offset_mm{};
    double angle_deg{8.0};
    double edge_fraction{0.82};
    double point_source_z_mm{-45.0};
    double point_source_y_mm{};
};

struct EditorState {
    std::vector<SurfaceEdit> surfaces;
    RayControls rays;
    osys::OpticalPresetId selected_preset{osys::OpticalPresetId::double_gauss_50mm_f2};
    bool show_layers_panel{true};
    bool show_rays_panel{true};
    bool show_presets_panel{true};
    bool show_view_panel{false};
    bool show_grid{true};
    bool show_aperture{true};
    bool show_mouse_debug{false};
    double zoom{1.0};
    double pan_z_mm{};
    double pan_y_mm{};
    int active_control{};
};

struct SourceRay {
    osys::Ray ray;
};

struct SpectralRay {
    osys::Ray ray;
    Color color{};
};

struct RayPath {
    std::vector<osys::Vec3> points;
    osys::TraceStatus status{osys::TraceStatus::ok};
};

struct WorldBounds {
    double min_z{-70.0};
    double max_z{90.0};
    double max_radius{30.0};
};

struct View2D {
    Rectangle canvas{};
    double center_z_mm{};
    double center_y_mm{};
    double px_per_mm{1.0};
};

struct UiFrame {
    Vector2 mouse{};
    int next_id{1};
    int* active_control{};
};

struct PanelCursor {
    float x{};
    float y{};
    float width{};
};

[[nodiscard]] double clamp_double(const double value, const double lo, const double hi) {
    return std::min(std::max(value, lo), hi);
}

[[nodiscard]] int clamp_int(const int value, const int lo, const int hi) {
    return std::min(std::max(value, lo), hi);
}

[[nodiscard]] double lerp_double(const double a, const double b, const double t) {
    return a + (b - a) * t;
}

[[nodiscard]] double deg_to_rad(const double degrees) {
    return degrees * k_pi / 180.0;
}

[[nodiscard]] bool near_zero(const double value) {
    return std::abs(value) < 1.0e-12;
}

[[nodiscard]] Vector2 ui_mouse_position() {
    return GetMousePosition();
}

[[nodiscard]] const char* surface_kind_name(const SurfaceKind kind) {
    switch (kind) {
        case SurfaceKind::plane:
            return "PlaneSagitta";
        case SurfaceKind::conic:
            return "ConicSagitta";
    }

    return "Unknown";
}

[[nodiscard]] const char* ray_preset_name(const RayPreset preset) {
    switch (preset) {
        case RayPreset::single_offset:
            return "Single offset";
        case RayPreset::parallel_bundle:
            return "Parallel bundle";
        case RayPreset::angled_center:
            return "Angled center";
        case RayPreset::edge_pair:
            return "Edge pair";
        case RayPreset::point_source:
            return "Point source";
    }

    return "Unknown";
}

[[nodiscard]] bool is_air_medium(const osys::Medium& medium) {
    return std::abs(osys::refractive_index(medium, 550.0) - 1.0) < 1.0e-6;
}

[[nodiscard]] Color medium_color(const osys::Medium& medium, const unsigned char alpha) {
    const double n = osys::refractive_index(medium, 550.0);
    if (n < 1.01) {
        return Color{120, 132, 146, alpha};
    }
    if (n < 1.56) {
        return Color{92, 176, 255, alpha};
    }
    if (n < 1.68) {
        return Color{96, 214, 188, alpha};
    }
    if (n < 1.78) {
        return Color{255, 183, 92, alpha};
    }

    return Color{226, 112, 168, alpha};
}

void set_surface_medium_from_catalog(SurfaceEdit& surface, const int catalog_index) {
    const int count = static_cast<int>(osys::medium_catalog_ids.size());
    const int index = (catalog_index % count + count) % count;
    const osys::MediumId id = osys::medium_catalog_ids[static_cast<std::size_t>(index)];
    surface.medium_after = osys::medium(id);
    surface.medium_name_after = std::string(osys::medium_name(id));
    surface.medium_catalog_index = index;
}

[[nodiscard]] osys::SagittaSurface to_sagitta_surface(const SurfaceEdit& edit) {
    if (edit.kind == SurfaceKind::plane) {
        return osys::SagittaSurface{osys::PlaneSagitta{}};
    }

    return osys::SagittaSurface{osys::ConicSagitta{
        .radius_mm = edit.radius_mm,
        .conic_constant = edit.conic_constant,
    }};
}

[[nodiscard]] osys::OpticalSurface to_optical_surface(const SurfaceEdit& edit) {
    return osys::OpticalSurface{
        .vertex_z_mm = edit.vertex_z_mm,
        .aperture_radius_mm = edit.aperture_radius_mm,
        .medium_after = edit.medium_after,
        .shape = to_sagitta_surface(edit),
    };
}

[[nodiscard]] osys::OpticalSystem to_optical_system(const std::vector<SurfaceEdit>& edits) {
    osys::OpticalSystem system{osys::air_medium};
    system.surfaces.reserve(edits.size());

    for (const SurfaceEdit& edit : edits) {
        osys::add_surface(system, to_optical_surface(edit));
    }

    return system;
}

[[nodiscard]] SurfaceEdit surface_edit_from_spec(const osys::OpticalSurfaceSpec& spec) {
    SurfaceEdit edit{
        .expanded = false,
        .vertex_z_mm = spec.vertex_z_mm,
        .aperture_radius_mm = spec.aperture_radius_mm,
        .medium_after = spec.medium_after,
        .medium_name_after = std::string(spec.medium_name_after),
        .medium_catalog_index = -1,
        .kind = SurfaceKind::plane,
        .radius_mm = 50.0,
        .conic_constant = 0.0,
    };

    if (std::holds_alternative<osys::ConicSagitta>(spec.shape.model)) {
        const osys::ConicSagitta conic = std::get<osys::ConicSagitta>(spec.shape.model);
        edit.kind = SurfaceKind::conic;
        edit.radius_mm = conic.radius_mm;
        edit.conic_constant = conic.conic_constant;
    }

    return edit;
}

[[nodiscard]] std::vector<SurfaceEdit> surfaces_for_preset(const osys::OpticalPresetId id) {
    const std::span<const osys::OpticalSurfaceSpec> specs = osys::optical_preset_surfaces(id);
    std::vector<SurfaceEdit> edits;
    edits.reserve(specs.size());
    for (const osys::OpticalSurfaceSpec& spec : specs) {
        edits.push_back(surface_edit_from_spec(spec));
    }
    return edits;
}

[[nodiscard]] std::vector<SurfaceEdit> default_surfaces() {
    return surfaces_for_preset(osys::OpticalPresetId::double_gauss_50mm_f2);
}

void load_optical_preset(EditorState& state, const osys::OpticalPresetId id) {
    state.selected_preset = id;
    state.surfaces = surfaces_for_preset(id);
    state.active_control = 0;
}

void cycle_optical_preset(EditorState& state, const bool forward) {
    load_optical_preset(
        state,
        forward
            ? osys::next_optical_preset_id(state.selected_preset)
            : osys::previous_optical_preset_id(state.selected_preset));
}

[[nodiscard]] double max_aperture_radius(const std::vector<SurfaceEdit>& surfaces) {
    double radius = 12.5;
    for (const SurfaceEdit& surface : surfaces) {
        radius = std::max(radius, surface.aperture_radius_mm);
    }
    return radius;
}

[[nodiscard]] double first_vertex_z(const std::vector<SurfaceEdit>& surfaces) {
    if (surfaces.empty()) {
        return 0.0;
    }

    return surfaces.front().vertex_z_mm;
}

[[nodiscard]] double last_vertex_z(const std::vector<SurfaceEdit>& surfaces) {
    if (surfaces.empty()) {
        return 20.0;
    }

    return surfaces.back().vertex_z_mm;
}

[[nodiscard]] WorldBounds world_bounds(const EditorState& state) {
    WorldBounds bounds{};
    const double aperture = max_aperture_radius(state.surfaces);
    const double first_z = first_vertex_z(state.surfaces);
    const double last_z = last_vertex_z(state.surfaces);

    bounds.min_z = first_z - 65.0;
    bounds.max_z = last_z + 95.0;
    bounds.max_radius = std::max(20.0, aperture * 1.7);

    if (state.rays.preset == RayPreset::point_source) {
        bounds.min_z = std::min(bounds.min_z, state.rays.point_source_z_mm - 10.0);
        bounds.max_radius = std::max(bounds.max_radius, std::abs(state.rays.point_source_y_mm) + aperture + 8.0);
    }

    if (state.rays.preset == RayPreset::angled_center) {
        const double distance = first_z - bounds.min_z;
        bounds.max_radius = std::max(bounds.max_radius, std::abs(std::tan(deg_to_rad(state.rays.angle_deg)) * distance) + 8.0);
    }

    return bounds;
}

[[nodiscard]] View2D make_view(const Rectangle canvas, const EditorState& state) {
    const WorldBounds bounds = world_bounds(state);
    const double z_span = std::max(20.0, bounds.max_z - bounds.min_z);
    const double y_span = std::max(10.0, 2.0 * bounds.max_radius);
    const double usable_w = std::max(100.0f, canvas.width - 80.0f);
    const double usable_h = std::max(100.0f, canvas.height - 80.0f);
    const double base_scale = std::min(usable_w / z_span, usable_h / y_span);

    return View2D{
        .canvas = canvas,
        .center_z_mm = 0.5 * (bounds.min_z + bounds.max_z) + state.pan_z_mm,
        .center_y_mm = state.pan_y_mm,
        .px_per_mm = base_scale * state.zoom,
    };
}

[[nodiscard]] Vector2 to_screen(const View2D& view, const osys::Vec3 point) {
    return Vector2{
        .x = view.canvas.x + 0.5f * view.canvas.width + static_cast<float>((point.z - view.center_z_mm) * view.px_per_mm),
        .y = view.canvas.y + 0.5f * view.canvas.height - static_cast<float>((point.y - view.center_y_mm) * view.px_per_mm),
    };
}

[[nodiscard]] osys::Vec3 normalize_or_forward(const osys::Vec3 value) {
    const double len = osys::length(value);
    if (near_zero(len) || !std::isfinite(len)) {
        return {0.0, 0.0, 1.0};
    }

    return value / len;
}

[[nodiscard]] Color wavelength_color(double wavelength_nm, const unsigned char alpha = 210) {
    wavelength_nm = clamp_double(wavelength_nm, 380.0, 700.0);

    double r = 0.0;
    double g = 0.0;
    double b = 0.0;

    if (wavelength_nm < 440.0) {
        r = -(wavelength_nm - 440.0) / (440.0 - 380.0);
        b = 1.0;
    } else if (wavelength_nm < 490.0) {
        g = (wavelength_nm - 440.0) / (490.0 - 440.0);
        b = 1.0;
    } else if (wavelength_nm < 510.0) {
        g = 1.0;
        b = -(wavelength_nm - 510.0) / (510.0 - 490.0);
    } else if (wavelength_nm < 580.0) {
        r = (wavelength_nm - 510.0) / (580.0 - 510.0);
        g = 1.0;
    } else if (wavelength_nm < 645.0) {
        r = 1.0;
        g = -(wavelength_nm - 645.0) / (645.0 - 580.0);
    } else {
        r = 1.0;
    }

    double factor = 1.0;
    if (wavelength_nm < 420.0) {
        factor = 0.3 + 0.7 * (wavelength_nm - 380.0) / (420.0 - 380.0);
    } else if (wavelength_nm > 645.0) {
        factor = 0.3 + 0.7 * (700.0 - wavelength_nm) / (700.0 - 645.0);
    }

    const double gamma = 0.8;
    const auto channel = [factor, gamma](const double value) -> unsigned char {
        if (value <= 0.0) {
            return 0;
        }

        return static_cast<unsigned char>(clamp_int(
            static_cast<int>(std::round(255.0 * std::pow(value * factor, gamma))),
            0,
            255));
    };

    return Color{channel(r), channel(g), channel(b), alpha};
}

[[nodiscard]] std::vector<double> wavelengths_for_controls(const RayControls& controls) {
    if (!controls.spectral) {
        return {controls.wavelength_nm};
    }

    const int sample_count = clamp_int(controls.spectral_samples, 3, 25);
    std::vector<double> wavelengths;
    wavelengths.reserve(static_cast<std::size_t>(sample_count));
    for (int i = 0; i < sample_count; ++i) {
        const double t = sample_count == 1 ? 0.0 : static_cast<double>(i) / static_cast<double>(sample_count - 1);
        wavelengths.push_back(lerp_double(450.0, 650.0, t));
    }
    return wavelengths;
}

[[nodiscard]] std::vector<SourceRay> source_rays_for_controls(const EditorState& state) {
    const double aperture = max_aperture_radius(state.surfaces);
    const double first_z = first_vertex_z(state.surfaces);
    const double start_z = first_z - 55.0;

    std::vector<SourceRay> source_rays;

    switch (state.rays.preset) {
        case RayPreset::single_offset: {
            source_rays.push_back(SourceRay{osys::Ray{
                .origin_mm = {0.0, state.rays.single_offset_mm, start_z},
                .direction = {0.0, 0.0, 1.0},
                .wavelength_nm = state.rays.wavelength_nm,
            }});
            break;
        }

        case RayPreset::parallel_bundle: {
            const int count = clamp_int(state.rays.ray_count, 1, 21);
            source_rays.reserve(static_cast<std::size_t>(count));
            for (int i = 0; i < count; ++i) {
                const double t = count == 1 ? 0.5 : static_cast<double>(i) / static_cast<double>(count - 1);
                const double y = lerp_double(-0.82 * aperture, 0.82 * aperture, t);
                source_rays.push_back(SourceRay{osys::Ray{
                    .origin_mm = {0.0, y, start_z},
                    .direction = {0.0, 0.0, 1.0},
                    .wavelength_nm = state.rays.wavelength_nm,
                }});
            }
            break;
        }

        case RayPreset::angled_center: {
            const double angle = deg_to_rad(state.rays.angle_deg);
            const double distance = first_z - start_z;
            const osys::Vec3 origin{0.0, -std::tan(angle) * distance, start_z};
            const osys::Vec3 target{0.0, 0.0, first_z};
            source_rays.push_back(SourceRay{osys::Ray{
                .origin_mm = origin,
                .direction = normalize_or_forward(target - origin),
                .wavelength_nm = state.rays.wavelength_nm,
            }});
            break;
        }

        case RayPreset::edge_pair: {
            const double y = clamp_double(state.rays.edge_fraction, 0.1, 0.98) * aperture;
            source_rays.push_back(SourceRay{osys::Ray{
                .origin_mm = {0.0, -y, start_z},
                .direction = {0.0, 0.0, 1.0},
                .wavelength_nm = state.rays.wavelength_nm,
            }});
            source_rays.push_back(SourceRay{osys::Ray{
                .origin_mm = {0.0, y, start_z},
                .direction = {0.0, 0.0, 1.0},
                .wavelength_nm = state.rays.wavelength_nm,
            }});
            break;
        }

        case RayPreset::point_source: {
            const int count = clamp_int(state.rays.point_samples, 2, 25);
            const osys::Vec3 origin{0.0, state.rays.point_source_y_mm, state.rays.point_source_z_mm};
            source_rays.reserve(static_cast<std::size_t>(count));
            for (int i = 0; i < count; ++i) {
                const double t = count == 1 ? 0.5 : static_cast<double>(i) / static_cast<double>(count - 1);
                const double target_y = lerp_double(-0.92 * aperture, 0.92 * aperture, t);
                const osys::Vec3 target{0.0, target_y, first_z};
                source_rays.push_back(SourceRay{osys::Ray{
                    .origin_mm = origin,
                    .direction = normalize_or_forward(target - origin),
                    .wavelength_nm = state.rays.wavelength_nm,
                }});
            }
            break;
        }
    }

    return source_rays;
}

[[nodiscard]] std::vector<SpectralRay> spectral_rays_for_controls(const EditorState& state) {
    const std::vector<SourceRay> sources = source_rays_for_controls(state);
    const std::vector<double> wavelengths = wavelengths_for_controls(state.rays);

    std::vector<SpectralRay> rays;
    rays.reserve(sources.size() * wavelengths.size());
    for (const SourceRay& source : sources) {
        for (const double wavelength : wavelengths) {
            osys::Ray ray = source.ray;
            ray.wavelength_nm = wavelength;
            rays.push_back(SpectralRay{
                .ray = ray,
                .color = wavelength_color(wavelength, state.rays.spectral ? 92 : 230),
            });
        }
    }

    return rays;
}

[[nodiscard]] std::optional<osys::Vec3> refract_direction(
    const osys::Vec3 incident,
    const osys::Vec3 surface_normal,
    const double n_incident,
    const double n_transmitted) {
    osys::Vec3 normal = osys::normalize(surface_normal);
    const osys::Vec3 direction = osys::normalize(incident);

    if (osys::dot(direction, normal) > 0.0) {
        normal = -normal;
    }

    const double eta = n_incident / n_transmitted;
    const double cos_i = -osys::dot(normal, direction);
    const double k = 1.0 - eta * eta * (1.0 - cos_i * cos_i);
    if (k < 0.0) {
        return std::nullopt;
    }

    return osys::normalize(eta * direction + (eta * cos_i - std::sqrt(k)) * normal);
}

[[nodiscard]] RayPath trace_path(
    const osys::OpticalSystem& system,
    const osys::Ray& input,
    const double end_z_mm) {
    RayPath path{};
    osys::Ray ray = input;
    ray.direction = normalize_or_forward(ray.direction);
    osys::Medium current_medium = system.initial_medium;
    path.points.push_back(ray.origin_mm);

    for (const osys::OpticalSurface& surface : system.surfaces) {
        const std::optional<double> hit_t = osys::intersect_surface(ray, surface);
        if (!hit_t.has_value()) {
            const double fallback_t = std::abs(end_z_mm - ray.origin_mm.z) + 40.0;
            path.points.push_back(ray.origin_mm + fallback_t * ray.direction);
            path.status = osys::TraceStatus::no_intersection;
            return path;
        }

        const osys::Vec3 hit = ray.origin_mm + *hit_t * ray.direction;
        path.points.push_back(hit);

        if (std::hypot(hit.x, hit.y) > surface.aperture_radius_mm) {
            path.status = osys::TraceStatus::missed_aperture;
            return path;
        }

        const std::optional<osys::Vec3> normal = osys::surface_normal(hit, surface);
        if (!normal.has_value()) {
            path.status = osys::TraceStatus::no_intersection;
            return path;
        }

        const double n_before = osys::refractive_index(current_medium, ray.wavelength_nm);
        const double n_after = osys::refractive_index(surface.medium_after, ray.wavelength_nm);
        const std::optional<osys::Vec3> refracted = refract_direction(ray.direction, *normal, n_before, n_after);
        if (!refracted.has_value()) {
            path.status = osys::TraceStatus::total_internal_reflection;
            return path;
        }

        ray.origin_mm = hit + k_ray_epsilon_mm * *refracted;
        ray.direction = *refracted;
        current_medium = surface.medium_after;
    }

    double t = 70.0;
    if (!near_zero(ray.direction.z)) {
        t = (end_z_mm - ray.origin_mm.z) / ray.direction.z;
        if (t < 1.0 || !std::isfinite(t)) {
            t = 70.0;
        }
    }

    path.points.push_back(ray.origin_mm + t * ray.direction);
    path.status = osys::TraceStatus::ok;
    return path;
}

[[nodiscard]] Rectangle row_rect(PanelCursor& cursor, const float height, const float gap = 6.0f) {
    Rectangle rect{
        .x = cursor.x,
        .y = cursor.y,
        .width = cursor.width,
        .height = height,
    };
    cursor.y += height + gap;
    return rect;
}

[[nodiscard]] bool point_in_rect(const Vector2 point, const Rectangle rect) {
    return CheckCollisionPointRec(point, rect);
}

[[nodiscard]] int next_control_id(UiFrame& ui) {
    const int id = ui.next_id;
    ++ui.next_id;
    return id;
}

void draw_panel_background(const Rectangle panel, const char* title) {
    DrawRectangleRounded(panel, 0.035f, 8, Color{22, 25, 31, 232});
    DrawRectangleRoundedLines(panel, 0.035f, 8, Color{88, 96, 112, 190});
    DrawText(title, static_cast<int>(panel.x + 14.0f), static_cast<int>(panel.y + 12.0f), 18, Color{235, 238, 244, 255});
}

bool ui_button(UiFrame& ui, const Rectangle rect, const char* label) {
    const bool hovered = point_in_rect(ui.mouse, rect);
    const bool down = hovered && IsMouseButtonDown(MOUSE_BUTTON_LEFT);
    const Color fill = down ? Color{80, 126, 180, 255}
        : hovered ? Color{62, 73, 92, 255}
        : Color{42, 48, 60, 255};

    DrawRectangleRounded(rect, 0.16f, 6, fill);
    DrawRectangleRoundedLines(rect, 0.16f, 6, Color{100, 112, 132, 180});

    const int font_size = 15;
    const int text_w = MeasureText(label, font_size);
    DrawText(
        label,
        static_cast<int>(rect.x + 0.5f * (rect.width - static_cast<float>(text_w))),
        static_cast<int>(rect.y + 0.5f * (rect.height - static_cast<float>(font_size))),
        font_size,
        Color{236, 239, 245, 255});

    return hovered && IsMouseButtonPressed(MOUSE_BUTTON_LEFT);
}

bool ui_toggle(UiFrame& ui, const Rectangle rect, const char* label, bool& value) {
    const bool clicked = ui_button(ui, rect, label);
    if (clicked) {
        value = !value;
    }

    if (value) {
        DrawRectangleRoundedLines(rect, 0.16f, 6, Color{120, 202, 255, 255});
    }

    return clicked;
}

bool ui_panel_close_button(UiFrame& ui, const Rectangle panel) {
    return ui_button(ui, Rectangle{panel.x + panel.width - 38.0f, panel.y + 8.0f, 24.0f, 24.0f}, "x");
}

bool ui_cycle_medium(UiFrame& ui, const Rectangle rect, SurfaceEdit& surface) {
    if (!ui_button(ui, rect, surface.medium_name_after.c_str())) {
        return false;
    }

    const int next = surface.medium_catalog_index < 0 ? 0 : surface.medium_catalog_index + 1;
    set_surface_medium_from_catalog(surface, next);
    return true;
}

bool ui_cycle_surface_kind(UiFrame& ui, const Rectangle rect, SurfaceKind& value) {
    if (!ui_button(ui, rect, surface_kind_name(value))) {
        return false;
    }

    value = value == SurfaceKind::plane ? SurfaceKind::conic : SurfaceKind::plane;
    return true;
}

bool ui_cycle_ray_preset(UiFrame& ui, const Rectangle rect, RayPreset& value) {
    if (!ui_button(ui, rect, ray_preset_name(value))) {
        return false;
    }

    switch (value) {
        case RayPreset::single_offset:
            value = RayPreset::parallel_bundle;
            break;
        case RayPreset::parallel_bundle:
            value = RayPreset::angled_center;
            break;
        case RayPreset::angled_center:
            value = RayPreset::edge_pair;
            break;
        case RayPreset::edge_pair:
            value = RayPreset::point_source;
            break;
        case RayPreset::point_source:
            value = RayPreset::single_offset;
            break;
    }

    return true;
}

void ui_label_value(const Rectangle rect, const char* label, const char* value) {
    DrawText(label, static_cast<int>(rect.x), static_cast<int>(rect.y + 2.0f), 14, Color{177, 185, 198, 255});
    DrawText(value, static_cast<int>(rect.x + rect.width - static_cast<float>(MeasureText(value, 14))), static_cast<int>(rect.y + 2.0f), 14, Color{232, 236, 242, 255});
}

bool ui_slider_double(
    UiFrame& ui,
    const Rectangle rect,
    const char* label,
    double& value,
    const double min_value,
    const double max_value,
    const char* format = "%.2f") {
    const int id = next_control_id(ui);
    const Rectangle label_rect{rect.x, rect.y, rect.width, 18.0f};
    const Rectangle track_rect{rect.x, rect.y + 22.0f, rect.width, 18.0f};
    char buffer[64]{};
    std::snprintf(buffer, sizeof(buffer), format, value);
    ui_label_value(label_rect, label, buffer);

    const bool hovered = point_in_rect(ui.mouse, track_rect);
    if (hovered && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
        *ui.active_control = id;
    }

    bool changed = false;
    if (*ui.active_control == id) {
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            const double t = clamp_double((ui.mouse.x - track_rect.x) / track_rect.width, 0.0, 1.0);
            value = lerp_double(min_value, max_value, t);
            changed = true;
        } else {
            *ui.active_control = 0;
        }
    }

    value = clamp_double(value, min_value, max_value);
    const float fill_w = static_cast<float>((value - min_value) / (max_value - min_value)) * track_rect.width;
    DrawRectangleRounded(track_rect, 0.5f, 8, Color{38, 43, 53, 255});
    DrawRectangleRounded(Rectangle{track_rect.x, track_rect.y, fill_w, track_rect.height}, 0.5f, 8, Color{77, 148, 219, 255});
    DrawCircleV(Vector2{track_rect.x + fill_w, track_rect.y + 0.5f * track_rect.height}, 7.0f, hovered || *ui.active_control == id ? Color{248, 250, 252, 255} : Color{195, 204, 216, 255});

    return changed;
}

bool ui_slider_int(
    UiFrame& ui,
    const Rectangle rect,
    const char* label,
    int& value,
    const int min_value,
    const int max_value) {
    double temp = static_cast<double>(value);
    const bool changed = ui_slider_double(ui, rect, label, temp, static_cast<double>(min_value), static_cast<double>(max_value), "%.0f");
    value = clamp_int(static_cast<int>(std::round(temp)), min_value, max_value);
    return changed;
}

void draw_grid(const View2D& view) {
    const double grid_mm = 10.0;
    const double left_z = view.center_z_mm - 0.5 * view.canvas.width / view.px_per_mm;
    const double right_z = view.center_z_mm + 0.5 * view.canvas.width / view.px_per_mm;
    const double bottom_y = view.center_y_mm - 0.5 * view.canvas.height / view.px_per_mm;
    const double top_y = view.center_y_mm + 0.5 * view.canvas.height / view.px_per_mm;

    const int z_start = static_cast<int>(std::floor(left_z / grid_mm)) - 1;
    const int z_end = static_cast<int>(std::ceil(right_z / grid_mm)) + 1;
    const int y_start = static_cast<int>(std::floor(bottom_y / grid_mm)) - 1;
    const int y_end = static_cast<int>(std::ceil(top_y / grid_mm)) + 1;

    for (int z_step = z_start; z_step <= z_end; ++z_step) {
        const double z = z_step * grid_mm;
        const Vector2 a = to_screen(view, osys::Vec3{0.0, bottom_y, z});
        const Vector2 b = to_screen(view, osys::Vec3{0.0, top_y, z});
        DrawLineV(a, b, z_step == 0 ? Color{95, 105, 122, 180} : Color{43, 49, 60, 155});
    }

    for (int y_step = y_start; y_step <= y_end; ++y_step) {
        const double y = y_step * grid_mm;
        const Vector2 a = to_screen(view, osys::Vec3{0.0, y, left_z});
        const Vector2 b = to_screen(view, osys::Vec3{0.0, y, right_z});
        DrawLineV(a, b, y_step == 0 ? Color{95, 105, 122, 180} : Color{43, 49, 60, 155});
    }
}

[[nodiscard]] osys::Vec3 surface_point_at_y(const SurfaceEdit& edit, const double y_mm) {
    const osys::SagittaSurface surface = to_sagitta_surface(edit);
    const std::optional<double> sagitta = osys::sagitta_mm(surface, std::abs(y_mm));
    return osys::Vec3{
        .x = 0.0,
        .y = y_mm,
        .z = edit.vertex_z_mm + sagitta.value_or(0.0),
    };
}

void draw_surface_curve(const View2D& view, const SurfaceEdit& edit, const Color color) {
    constexpr int samples = 72;
    const double aperture = edit.aperture_radius_mm;
    Vector2 previous{};
    bool has_previous = false;

    for (int i = 0; i <= samples; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(samples);
        const double y = lerp_double(-aperture, aperture, t);
        const osys::Vec3 point = surface_point_at_y(edit, y);
        const Vector2 screen = to_screen(view, point);
        if (has_previous) {
            DrawLineEx(previous, screen, 2.2f, color);
        }
        previous = screen;
        has_previous = true;
    }
}

void draw_lens_medium_regions(const View2D& view, const std::vector<SurfaceEdit>& surfaces) {
    if (surfaces.size() < 2) {
        return;
    }

    constexpr int samples = 48;
    for (std::size_t i = 0; i + 1 < surfaces.size(); ++i) {
        const SurfaceEdit& left = surfaces[i];
        const SurfaceEdit& right = surfaces[i + 1];
        if (is_air_medium(left.medium_after)) {
            continue;
        }

        const double aperture = std::min(left.aperture_radius_mm, right.aperture_radius_mm);
        const Color fill = medium_color(left.medium_after, 38);

        for (int sample = 0; sample < samples; ++sample) {
            const double t0 = static_cast<double>(sample) / static_cast<double>(samples);
            const double t1 = static_cast<double>(sample + 1) / static_cast<double>(samples);
            const double y0 = lerp_double(-aperture, aperture, t0);
            const double y1 = lerp_double(-aperture, aperture, t1);

            const Vector2 l0 = to_screen(view, surface_point_at_y(left, y0));
            const Vector2 r0 = to_screen(view, surface_point_at_y(right, y0));
            const Vector2 l1 = to_screen(view, surface_point_at_y(left, y1));
            const Vector2 r1 = to_screen(view, surface_point_at_y(right, y1));

            DrawTriangle(l0, r0, l1, fill);
            DrawTriangle(l1, r0, r1, fill);
        }
    }
}

void draw_aperture_lines(const View2D& view, const std::vector<SurfaceEdit>& surfaces) {
    for (const SurfaceEdit& surface : surfaces) {
        const Vector2 top = to_screen(view, osys::Vec3{0.0, surface.aperture_radius_mm, surface.vertex_z_mm});
        const Vector2 bottom = to_screen(view, osys::Vec3{0.0, -surface.aperture_radius_mm, surface.vertex_z_mm});
        DrawCircleV(top, 3.0f, Color{214, 221, 232, 180});
        DrawCircleV(bottom, 3.0f, Color{214, 221, 232, 180});
    }
}

void draw_surfaces(const View2D& view, const std::vector<SurfaceEdit>& surfaces, const bool show_aperture) {
    draw_lens_medium_regions(view, surfaces);

    for (const SurfaceEdit& surface : surfaces) {
        draw_surface_curve(view, surface, medium_color(surface.medium_after, 235));
    }

    if (show_aperture) {
        draw_aperture_lines(view, surfaces);
    }
}

void draw_ray_path(const View2D& view, const RayPath& path, const Color color, const float width) {
    if (path.points.size() < 2) {
        return;
    }

    for (std::size_t i = 1; i < path.points.size(); ++i) {
        const Vector2 a = to_screen(view, path.points[i - 1]);
        const Vector2 b = to_screen(view, path.points[i]);
        DrawLineEx(a, b, width, color);
    }

    if (path.status != osys::TraceStatus::ok) {
        const Vector2 last = to_screen(view, path.points.back());
        DrawCircleLines(static_cast<int>(last.x), static_cast<int>(last.y), 5.0f, Color{255, 96, 96, 230});
    }
}

void draw_rays(const View2D& view, const EditorState& state) {
    const osys::OpticalSystem system = to_optical_system(state.surfaces);
    const std::vector<SpectralRay> rays = spectral_rays_for_controls(state);
    const double end_z = world_bounds(state).max_z;

    BeginBlendMode(state.rays.spectral ? BLEND_ADDITIVE : BLEND_ALPHA);
    for (const SpectralRay& ray : rays) {
        const RayPath path = trace_path(system, ray.ray, end_z);
        draw_ray_path(view, path, ray.color, state.rays.spectral ? 1.2f : 1.8f);
    }
    EndBlendMode();
}

void draw_canvas(const EditorState& state) {
    const Rectangle canvas{0.0f, 0.0f, static_cast<float>(GetScreenWidth()), static_cast<float>(GetScreenHeight())};
    const View2D view = make_view(canvas, state);

    ClearBackground(Color{9, 11, 16, 255});
    if (state.show_grid) {
        draw_grid(view);
    }

    draw_rays(view, state);
    draw_surfaces(view, state.surfaces, state.show_aperture);
}

void add_surface_after_last(std::vector<SurfaceEdit>& surfaces) {
    const double next_z = surfaces.empty() ? 0.0 : surfaces.back().vertex_z_mm + 5.0;
    const double aperture = surfaces.empty() ? 12.5 : surfaces.back().aperture_radius_mm;
    const double radius = surfaces.empty() ? 50.0 : -surfaces.back().radius_mm;
    SurfaceEdit surface{
        .expanded = true,
        .vertex_z_mm = next_z,
        .aperture_radius_mm = aperture,
        .kind = SurfaceKind::conic,
        .radius_mm = near_zero(radius) ? 50.0 : radius,
        .conic_constant = 0.0,
    };
    set_surface_medium_from_catalog(surface, surfaces.empty() || !is_air_medium(surfaces.back().medium_after) ? 0 : 1);
    surfaces.push_back(surface);
}

void draw_layers_panel(UiFrame& ui, EditorState& state) {
    const float panel_w = 360.0f;
    const Rectangle panel{
        static_cast<float>(GetScreenWidth()) - panel_w - 18.0f,
        58.0f,
        panel_w,
        static_cast<float>(GetScreenHeight()) - 78.0f,
    };

    draw_panel_background(panel, "Layers");
    if (ui_panel_close_button(ui, panel)) {
        state.show_layers_panel = false;
        return;
    }

    PanelCursor cursor{panel.x + 14.0f, panel.y + 44.0f, panel.width - 28.0f};

    Rectangle header = row_rect(cursor, 30.0f);
    if (ui_button(ui, Rectangle{header.x, header.y, 46.0f, header.height}, "+")) {
        add_surface_after_last(state.surfaces);
    }
    if (ui_button(ui, Rectangle{header.x + 54.0f, header.y, 46.0f, header.height}, "-") && !state.surfaces.empty()) {
        state.surfaces.pop_back();
    }
    ui_toggle(ui, Rectangle{header.x + header.width - 112.0f, header.y, 112.0f, header.height}, "Aperture", state.show_aperture);

    for (std::size_t i = 0; i < state.surfaces.size(); ++i) {
        SurfaceEdit& surface = state.surfaces[i];
        const Rectangle row = row_rect(cursor, 28.0f, 3.0f);
        DrawRectangleRounded(row, 0.08f, 5, Color{34, 39, 49, 240});

        const std::string title = std::string(surface.expanded ? "v " : "> ") + "Surface " + std::to_string(i);
        if (ui_button(ui, Rectangle{row.x + 6.0f, row.y + 3.0f, 138.0f, 22.0f}, title.c_str())) {
            surface.expanded = !surface.expanded;
        }

        const std::string z_text = TextFormat("z %.2f", surface.vertex_z_mm);
        DrawText(z_text.c_str(), static_cast<int>(row.x + 154.0f), static_cast<int>(row.y + 7.0f), 13, Color{190, 199, 212, 255});

        if (ui_button(ui, Rectangle{row.x + row.width - 34.0f, row.y + 3.0f, 28.0f, 22.0f}, "x")) {
            state.surfaces.erase(state.surfaces.begin() + static_cast<std::ptrdiff_t>(i));
            break;
        }

        if (!surface.expanded) {
            continue;
        }

        const Rectangle type_row = row_rect(cursor, 28.0f);
        DrawText("Type", static_cast<int>(type_row.x), static_cast<int>(type_row.y + 7.0f), 14, Color{177, 185, 198, 255});
        ui_cycle_surface_kind(ui, Rectangle{type_row.x + 112.0f, type_row.y, type_row.width - 112.0f, type_row.height}, surface.kind);

        const Rectangle medium_row = row_rect(cursor, 28.0f);
        DrawText("Medium after", static_cast<int>(medium_row.x), static_cast<int>(medium_row.y + 7.0f), 14, Color{177, 185, 198, 255});
        ui_cycle_medium(ui, Rectangle{medium_row.x + 112.0f, medium_row.y, medium_row.width - 112.0f, medium_row.height}, surface);

        ui_slider_double(ui, row_rect(cursor, 44.0f), "vertex z mm", surface.vertex_z_mm, -60.0, 120.0);
        ui_slider_double(ui, row_rect(cursor, 44.0f), "aperture radius mm", surface.aperture_radius_mm, 1.0, 35.0);

        if (surface.kind == SurfaceKind::conic) {
            const double radius_limit = std::max(160.0, std::ceil(std::abs(surface.radius_mm) / 50.0) * 50.0);
            ui_slider_double(ui, row_rect(cursor, 44.0f), "radius mm", surface.radius_mm, -radius_limit, radius_limit);
            if (std::abs(surface.radius_mm) < 1.0) {
                surface.radius_mm = surface.radius_mm < 0.0 ? -1.0 : 1.0;
            }
            ui_slider_double(ui, row_rect(cursor, 44.0f), "conic constant", surface.conic_constant, -2.5, 2.5);
        }

        cursor.y += 8.0f;
    }
}

void draw_rays_panel(UiFrame& ui, EditorState& state) {
    const float panel_w = 360.0f;
    const float x = state.show_layers_panel
        ? static_cast<float>(GetScreenWidth()) - 2.0f * panel_w - 30.0f
        : static_cast<float>(GetScreenWidth()) - panel_w - 18.0f;
    const Rectangle panel{x, 58.0f, panel_w, 536.0f};

    draw_panel_background(panel, "Rays");
    if (ui_panel_close_button(ui, panel)) {
        state.show_rays_panel = false;
        return;
    }

    PanelCursor cursor{panel.x + 14.0f, panel.y + 44.0f, panel.width - 28.0f};

    const Rectangle preset_row = row_rect(cursor, 30.0f);
    DrawText("Preset", static_cast<int>(preset_row.x), static_cast<int>(preset_row.y + 8.0f), 14, Color{177, 185, 198, 255});
    ui_cycle_ray_preset(ui, Rectangle{preset_row.x + 106.0f, preset_row.y, preset_row.width - 106.0f, preset_row.height}, state.rays.preset);

    const Rectangle spectral_row = row_rect(cursor, 30.0f);
    ui_toggle(ui, Rectangle{spectral_row.x, spectral_row.y, 118.0f, spectral_row.height}, "Spectral", state.rays.spectral);
    if (!state.rays.spectral) {
        ui_slider_double(ui, row_rect(cursor, 44.0f), "wavelength nm", state.rays.wavelength_nm, 450.0, 650.0, "%.0f");
    } else {
        ui_slider_int(ui, row_rect(cursor, 44.0f), "spectral samples", state.rays.spectral_samples, 3, 25);
    }

    switch (state.rays.preset) {
        case RayPreset::single_offset:
            ui_slider_double(ui, row_rect(cursor, 44.0f), "ray offset mm", state.rays.single_offset_mm, -24.0, 24.0);
            break;
        case RayPreset::parallel_bundle:
            ui_slider_int(ui, row_rect(cursor, 44.0f), "input lines", state.rays.ray_count, 1, 21);
            break;
        case RayPreset::angled_center:
            ui_slider_double(ui, row_rect(cursor, 44.0f), "angle degrees", state.rays.angle_deg, -35.0, 35.0);
            break;
        case RayPreset::edge_pair:
            ui_slider_double(ui, row_rect(cursor, 44.0f), "edge fraction", state.rays.edge_fraction, 0.1, 0.98);
            break;
        case RayPreset::point_source:
            ui_slider_int(ui, row_rect(cursor, 44.0f), "surface samples", state.rays.point_samples, 2, 25);
            ui_slider_double(ui, row_rect(cursor, 44.0f), "source z mm", state.rays.point_source_z_mm, -110.0, -5.0);
            ui_slider_double(ui, row_rect(cursor, 44.0f), "source y mm", state.rays.point_source_y_mm, -35.0, 35.0);
            break;
    }
}

void draw_presets_panel(UiFrame& ui, EditorState& state) {
    const Rectangle panel{18.0f, 58.0f, 390.0f, 158.0f};
    draw_panel_background(panel, "Presets");
    if (ui_panel_close_button(ui, panel)) {
        state.show_presets_panel = false;
        return;
    }

    const osys::OpticalPresetInfo& info = osys::optical_preset_info(state.selected_preset);
    PanelCursor cursor{panel.x + 14.0f, panel.y + 44.0f, panel.width - 28.0f};

    const Rectangle name_row = row_rect(cursor, 34.0f);
    DrawText(info.name.data(), static_cast<int>(name_row.x), static_cast<int>(name_row.y + 8.0f), 15, Color{232, 236, 242, 255});

    const Rectangle meta_row = row_rect(cursor, 20.0f);
    const std::string meta = TextFormat("%.2f mm  f/%.2f  %zu surfaces", info.focal_length_mm, info.f_number, state.surfaces.size());
    DrawText(meta.c_str(), static_cast<int>(meta_row.x), static_cast<int>(meta_row.y + 2.0f), 14, Color{177, 185, 198, 255});

    const Rectangle action_row = row_rect(cursor, 30.0f);
    if (ui_button(ui, Rectangle{action_row.x, action_row.y, 84.0f, action_row.height}, "Prev")) {
        cycle_optical_preset(state, false);
    }
    if (ui_button(ui, Rectangle{action_row.x + 92.0f, action_row.y, 84.0f, action_row.height}, "Next")) {
        cycle_optical_preset(state, true);
    }
    if (ui_button(ui, Rectangle{action_row.x + 184.0f, action_row.y, 94.0f, action_row.height}, "Reload")) {
        load_optical_preset(state, state.selected_preset);
    }
}

void draw_view_panel(UiFrame& ui, EditorState& state) {
    const float panel_y = state.show_presets_panel ? 232.0f : 58.0f;
    const Rectangle panel{18.0f, panel_y, 306.0f, 232.0f};
    draw_panel_background(panel, "View");
    if (ui_panel_close_button(ui, panel)) {
        state.show_view_panel = false;
        return;
    }

    PanelCursor cursor{panel.x + 14.0f, panel.y + 44.0f, panel.width - 28.0f};

    ui_toggle(ui, row_rect(cursor, 30.0f), "Grid", state.show_grid);
    ui_slider_double(ui, row_rect(cursor, 44.0f), "zoom", state.zoom, 0.35, 2.8);
    ui_slider_double(ui, row_rect(cursor, 44.0f), "pan z mm", state.pan_z_mm, -80.0, 80.0);
    ui_slider_double(ui, row_rect(cursor, 44.0f), "pan y mm", state.pan_y_mm, -40.0, 40.0);

    const Rectangle reset_row = row_rect(cursor, 30.0f);
    if (ui_button(ui, reset_row, "Reset view")) {
        state.zoom = 1.0;
        state.pan_z_mm = 0.0;
        state.pan_y_mm = 0.0;
    }
}

void draw_toolbar(UiFrame& ui, EditorState& state) {
    const Rectangle bar{18.0f, 14.0f, 572.0f, 34.0f};
    DrawRectangleRounded(bar, 0.25f, 8, Color{22, 25, 31, 232});
    DrawRectangleRoundedLines(bar, 0.25f, 8, Color{88, 96, 112, 160});

    ui_toggle(ui, Rectangle{bar.x + 8.0f, bar.y + 5.0f, 94.0f, 24.0f}, "Layers", state.show_layers_panel);
    ui_toggle(ui, Rectangle{bar.x + 110.0f, bar.y + 5.0f, 78.0f, 24.0f}, "Rays", state.show_rays_panel);
    ui_toggle(ui, Rectangle{bar.x + 196.0f, bar.y + 5.0f, 94.0f, 24.0f}, "Presets", state.show_presets_panel);
    ui_toggle(ui, Rectangle{bar.x + 298.0f, bar.y + 5.0f, 78.0f, 24.0f}, "View", state.show_view_panel);
    ui_toggle(ui, Rectangle{bar.x + 384.0f, bar.y + 5.0f, 78.0f, 24.0f}, "Grid", state.show_grid);

    if (ui_button(ui, Rectangle{bar.x + 470.0f, bar.y + 5.0f, 94.0f, 24.0f}, "Reload")) {
        load_optical_preset(state, state.selected_preset);
    }
}

void draw_status_footer(const EditorState& state) {
    const std::vector<SpectralRay> rays = spectral_rays_for_controls(state);
    const osys::OpticalPresetInfo& info = osys::optical_preset_info(state.selected_preset);
    const std::string text = TextFormat(
        "%s  |  %zu surfaces  |  %zu traced rays  |  %s",
        info.name.data(),
        state.surfaces.size(),
        rays.size(),
        state.rays.spectral ? "additive spectral blend" : "single wavelength");
    DrawRectangle(18, GetScreenHeight() - 38, MeasureText(text.c_str(), 14) + 24, 24, Color{12, 15, 22, 190});
    DrawText(text.c_str(), 30, GetScreenHeight() - 32, 14, Color{190, 199, 212, 255});
}

void draw_mouse_debug(const Vector2 ui_mouse) {
    const Vector2 raw = GetMousePosition();
    const Vector2 dpi = GetWindowScaleDPI();
    const std::string text = TextFormat(
        "mouse raw %.0f,%.0f  ui %.0f,%.0f  screen %dx%d  render %dx%d  dpi %.2f,%.2f",
        raw.x,
        raw.y,
        ui_mouse.x,
        ui_mouse.y,
        GetScreenWidth(),
        GetScreenHeight(),
        GetRenderWidth(),
        GetRenderHeight(),
        dpi.x,
        dpi.y);

    DrawLine(static_cast<int>(ui_mouse.x - 8.0f), static_cast<int>(ui_mouse.y), static_cast<int>(ui_mouse.x + 8.0f), static_cast<int>(ui_mouse.y), Color{255, 80, 80, 255});
    DrawLine(static_cast<int>(ui_mouse.x), static_cast<int>(ui_mouse.y - 8.0f), static_cast<int>(ui_mouse.x), static_cast<int>(ui_mouse.y + 8.0f), Color{255, 80, 80, 255});
    DrawRectangle(18, 296, MeasureText(text.c_str(), 14) + 24, 24, Color{12, 15, 22, 210});
    DrawText(text.c_str(), 30, 302, 14, Color{245, 222, 160, 255});
}

void draw_editor(EditorState& state) {
    if (IsKeyPressed(KEY_F1)) {
        state.show_mouse_debug = !state.show_mouse_debug;
    }

    draw_canvas(state);

    const Vector2 ui_mouse = ui_mouse_position();
    UiFrame ui{
        .mouse = ui_mouse,
        .next_id = 1,
        .active_control = &state.active_control,
    };

    draw_toolbar(ui, state);

    if (state.show_presets_panel) {
        draw_presets_panel(ui, state);
    }

    if (state.show_layers_panel) {
        draw_layers_panel(ui, state);
    }

    if (state.show_rays_panel) {
        draw_rays_panel(ui, state);
    }

    if (state.show_view_panel) {
        draw_view_panel(ui, state);
    }

    draw_status_footer(state);

    if (state.show_mouse_debug) {
        draw_mouse_debug(ui_mouse);
    }
}

} // namespace

int main() {
    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE | FLAG_VSYNC_HINT);
    InitWindow(k_window_width, k_window_height, "osys editor");
    SetWindowMinSize(k_min_window_width, k_min_window_height);
    SetTargetFPS(60);

    EditorState state{
        .surfaces = default_surfaces(),
    };

    while (!WindowShouldClose()) {
        BeginDrawing();
        draw_editor(state);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
