#pragma once

#include <Eigen/Core>
#include <optional>

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

    bool is_nearly_quadratic() const;
    bool is_nearly_linear() const;
    bool is_nearly_constant() const;
};

struct RootInterval {
    RootInterval(const double a, const double b, const double min_gradient)
        : a(a)
        , b(b)
        , min_gradient(min_gradient)
    {
    }

    double a, b;
    double min_gradient;
};

/// @brief Compute the time of impact between a point and triangle.
/// @param p_t0 Point at time 0.
/// @param t0_t0 Triangle vertex 0 at time 0.
/// @param t1_t0 Triangle vertex 1 at time 0.
/// @param t2_t0 Triangle vertex 2 at time 0.
/// @param p_t1 Point at time 1.
/// @param t0_t1 Triangle vertex 0 at time 1.
/// @param t1_t1 Triangle vertex 1 at time 1.
/// @param t2_t1 Triangle vertex 2 at time 1.
/// @param toi Computed time of impact.
/// @return True if the point and triangle collide, false otherwise.
bool point_triangle_ccd(
    const Eigen::Ref<const Eigen::Vector3d>& p_t0,
    const Eigen::Ref<const Eigen::Vector3d>& t0_t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1_t0,
    const Eigen::Ref<const Eigen::Vector3d>& t2_t0,
    const Eigen::Ref<const Eigen::Vector3d>& p_t1,
    const Eigen::Ref<const Eigen::Vector3d>& t0_t1,
    const Eigen::Ref<const Eigen::Vector3d>& t1_t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2_t1,
    double& toi);

/// @brief Compute the time of impact between two edges.
/// @param ea0_t0 Edge A vertex 0 at time 0.
/// @param ea1_t0 Edge A vertex 1 at time 0.
/// @param eb0_t0 Edge B vertex 0 at time 0.
/// @param eb1_t0 Edge B vertex 1 at time 0.
/// @param ea0_t1 Edge A vertex 0 at time 1.
/// @param ea1_t1 Edge A vertex 1 at time 1.
/// @param eb0_t1 Edge B vertex 0 at time 1.
/// @param eb1_t1 Edge B vertex 1 at time 1.
/// @param toi Computed time of impact.
/// @return True if the edges collide, false otherwise.
bool edge_edge_ccd(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_t0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_t0,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_t0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_t0,
    const Eigen::Ref<const Eigen::Vector3d>& ea0_t1,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_t1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_t1,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_t1,
    double& toi);

/// @brief Determine the root interval of a cubic equation.
/// @param d Cubic equation to determine the root interval of.
/// @return Root interval of the cubic equation.
std::optional<RootInterval>
determine_cubic_root_interval(const CubicEquation& d);

/// @brief Perform a modified Newton-Raphson root finding algorithm to find the roots of a cubic equation.
/// @param f Cubic equation to find the roots of.
/// @param x0 Initial guess for the root.
/// @param tolerance Tolerance for the root finding algorithm.
/// @param max_iterations Maximum number of iterations for the root finding algorithm.
/// @return Root of the cubic equation.
double modified_newton_raphson(
    const CubicEquation& f,
    const double x0,
    const double locally_min_gradient,
    const double tolerance = 1e-8,
    const unsigned max_iterations = std::numeric_limits<unsigned>::max());

} // namespace ccd