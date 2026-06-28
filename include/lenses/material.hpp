#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>

namespace lenses {

struct ConstantIndex {
    double n{};
};

struct Sellmeier3 {
    // C coefficients are in micrometer^2. Input wavelengths are still nm.
    std::array<double, 3> b{};
    std::array<double, 3> c_um2{};
};

struct Material {
    Material() = default;
    explicit Material(ConstantIndex constant) : model(constant) {}
    explicit Material(Sellmeier3 sellmeier) : model(sellmeier) {}

    [[nodiscard]] double refractive_index(double wavelength_nm) const;

    std::variant<ConstantIndex, Sellmeier3> model{ConstantIndex{1.0}};
};

[[nodiscard]] inline double refractive_index(const Material& material, double wavelength_nm) {
    return std::visit([wavelength_nm](const auto& model) {
        using Model = std::decay_t<decltype(model)>;

        if constexpr (std::is_same_v<Model, ConstantIndex>) {
            return model.n;
        } else {
            const double lambda_um = wavelength_nm * 0.001;
            const double lambda2 = lambda_um * lambda_um;
            double n2 = 1.0;

            for (std::size_t i = 0; i < model.b.size(); ++i) {
                n2 += model.b[i] * lambda2 / (lambda2 - model.c_um2[i]);
            }

            return std::sqrt(n2);
        }
    }, material.model);
}

inline double Material::refractive_index(double wavelength_nm) const {
    return lenses::refractive_index(*this, wavelength_nm);
}

struct MaterialCatalog {
    void add(std::string id, Material material);

    [[nodiscard]] const Material* find(const std::string& id) const;

    [[nodiscard]] double refractive_index(const std::string& id, double wavelength_nm) const;

    std::unordered_map<std::string, Material> materials;
};

inline void add_material(MaterialCatalog& catalog, std::string id, Material material) {
    catalog.materials.insert_or_assign(std::move(id), std::move(material));
}

[[nodiscard]] inline const Material* find_material(const MaterialCatalog& catalog, const std::string& id) {
    const auto found = catalog.materials.find(id);
    if (found == catalog.materials.end()) {
        return nullptr;
    }
    return &found->second;
}

[[nodiscard]] inline double refractive_index(
    const MaterialCatalog& catalog,
    const std::string& id,
    double wavelength_nm) {
    const Material* material = find_material(catalog, id);
    if (material == nullptr) {
        throw std::out_of_range("unknown optical material: " + id);
    }
    return refractive_index(*material, wavelength_nm);
}

inline void MaterialCatalog::add(std::string id, Material material) {
    add_material(*this, std::move(id), std::move(material));
}

inline const Material* MaterialCatalog::find(const std::string& id) const {
    return find_material(*this, id);
}

inline double MaterialCatalog::refractive_index(const std::string& id, double wavelength_nm) const {
    return lenses::refractive_index(*this, id, wavelength_nm);
}

[[nodiscard]] inline Material make_air() {
    return Material{ConstantIndex{1.0}};
}

[[nodiscard]] inline Material make_n_bk7() {
    return Material{Sellmeier3{
        .b = {1.03961212, 0.231792344, 1.01046945},
        .c_um2 = {0.00600069867, 0.0200179144, 103.560653},
    }};
}

} // namespace lenses
