#pragma once

#include "cubic.hpp"

#include <Eigen/Core>
#include <optional>

#include <functional>

namespace ccd {

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

/// @brief Compute the roots of a cubic equation.
/// @param d Cubic equation to compute the roots of.
/// @param toi Computed time of impact.
/// @return True if the cubic equation has a root in [0, 1], false otherwise.
bool fast_approximate_root_ccd(
    CubicEquation d,
    const std::function<bool(const double)>& is_inside,
    double& toi);

/// @brief Root interval of a polynomial.
struct RootInterval {
    /// @brief Construct a root interval.
    /// @param a Lower bound of the root interval.
    /// @param b Upper bound of the root interval.
    /// @param min_gradient  Minimum gradient of the polynomial in the interval.
    RootInterval(const double a, const double b, const double min_gradient)
        : a(a)
        , b(b)
        , min_gradient(min_gradient)
    {
    }

    double a; ///< Lower bound of the root interval.
    double b; ///< Upper bound of the root interval.
    /// Minimum gradient of the polynomial in the interval.
    double min_gradient;
};

/// @brief Determine the root interval of a cubic equation.
/// @param d Cubic equation to determine the root interval of.
/// @return Root interval of the cubic equation.
std::optional<RootInterval>
determine_cubic_root_interval(const CubicEquation& d);

} // namespace ccd