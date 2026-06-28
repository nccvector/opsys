#pragma once

#include <cmath>

namespace osys {

struct Vec3 {
    double x{};
    double y{};
    double z{};

    [[nodiscard]] constexpr Vec3 operator+(const Vec3 rhs) const {
        return {x + rhs.x, y + rhs.y, z + rhs.z};
    }

    [[nodiscard]] constexpr Vec3 operator-(const Vec3 rhs) const {
        return {x - rhs.x, y - rhs.y, z - rhs.z};
    }

    [[nodiscard]] constexpr Vec3 operator-() const {
        return {-x, -y, -z};
    }

    [[nodiscard]] constexpr Vec3 operator*(double scale) const {
        return {x * scale, y * scale, z * scale};
    }

    [[nodiscard]] constexpr Vec3 operator/(const double scale) const {
        return {x / scale, y / scale, z / scale};
    }
};

[[nodiscard]] constexpr Vec3 operator*(const double scale, const Vec3 &value) {
    return value * scale;
}

[[nodiscard]] constexpr double dot(const Vec3 &lhs, const Vec3 &rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

[[nodiscard]] inline double length(const Vec3 &value) {
    return std::sqrt(dot(value, value));
}

[[nodiscard]] inline Vec3 normalize(const Vec3 &value) {
    const double len = length(value);
    return value / len;
}

} // namespace osys
