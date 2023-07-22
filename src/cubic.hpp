#pragma once

#include <array>
#include <vector>
#include <limits>
#include <iostream>

namespace ccd {

/// @brief Cubic equation of the form ax³ + bx² + cx + d.
struct CubicEquation {
    /// @brief Coefficients of the cubic equation.
    double a, b, c, d;

    /// @brief Evaluate the cubic equation at t.
    /// @param t Value of t to evaluate the cubic equation at.
    /// @return Value of the cubic equation at t.
    double operator()(const double x) const
    {
        return x * (x * (x * a + b) + c) + d;
    }

    /// @brief Evaluate the derivative of the cubic equation at t.
    /// @param t Value of t to evaluate the derivative of the cubic equation at.
    /// @return Value of the derivative of the cubic equation at t.
    double derivative(const double x) const
    {
        return x * (x * 3 * a + 2 * b) + c;
    }

    std::array<double, 2> extrema() const;
    double inflection() const { return -b / (3 * a); }

    std::vector<std::array<double, 2>> monotonic_intervals() const;

    bool is_nearly_quadratic(
        const double tol = 10 * std::numeric_limits<double>::epsilon()) const;
    bool is_nearly_linear(
        const double tol = 10 * std::numeric_limits<double>::epsilon()) const;
    bool is_nearly_constant(
        const double tol = 10 * std::numeric_limits<double>::epsilon()) const;

    CubicEquation& operator*=(const double x)
    {
        this->a *= x;
        this->b *= x;
        this->c *= x;
        this->d *= x;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const CubicEquation& f);
};

} // namespace ccd