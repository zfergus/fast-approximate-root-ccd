#include "ccd.hpp"
#include "autogen.hpp"

namespace ccd {

bool point_triangle_ccd(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& t0_t0,
    const Eigen::Vector3d& t1_t0,
    const Eigen::Vector3d& t2_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& t0_t1,
    const Eigen::Vector3d& t1_t1,
    const Eigen::Vector3d& t2_t1,
    double& toi)
{
    const CubicEquation d = autogen::point_triangle_ccd_equation(
        p_t0.x(), p_t0.y(), p_t0.z(), t0_t0.x(), t0_t0.y(), t0_t0.z(),
        t1_t0.x(), t1_t0.y(), t1_t0.z(), t2_t0.x(), t2_t0.y(), t2_t0.z(),
        p_t1.x(), p_t1.y(), p_t1.z(), t0_t1.x(), t0_t1.y(), t0_t1.z(),
        t1_t1.x(), t1_t1.y(), t1_t1.z(), t2_t1.x(), t2_t1.y(), t2_t1.z());

    const auto root_interval = determine_cubic_root_interval(d);
    if (!root_interval || root_interval->a > root_interval->b) {
        return false;
    }

    toi = modified_newton_raphson(
        d, root_interval->a, root_interval->min_gradient);

    return toi >= 0 && toi <= 1;
}

std::optional<RootInterval>
determine_cubic_root_interval(const CubicEquation& d)
{
    double t0 = 0, t1 = 1;
    assert(d(t0) > 0);

    const auto [tm0, tm1] = d.extrema();
    assert(tm0 < tm1);

    assert(d.a != 0);
    if (d.a > 0) {
        // Case 1 (Interval II)
        const double t_min = tm1, t_max = tm0;
        if (d(t_min) > 0 || t0 >= t_min)
            return std::nullopt;
        return RootInterval(t_min, t_max, d.derivative(-d.a / (3 * d.a)));
    }

    const double t_min = tm0, t_max = tm1;
    if (d(t_min) >= 0) {
        // Case 2 (Interval III or may be an invalid interval if t_max > t1)
        return RootInterval(t_max, t1, d.derivative(t1));
    } else {
        // Case 3
        return (t0 < t_min)
            ? RootInterval(t0, t_min, d.derivative(t0))  // Interval I
            : RootInterval(t_max, t1, d.derivative(t1)); // Interval III
    }
}

double modified_newton_raphson(
    const CubicEquation& f,
    const double x0,
    const double locally_min_gradient,
    const double tolerance,
    const unsigned max_iterations)
{
    double xi = x0;
    double yi = f(xi);
    for (int i = 0; i < max_iterations && std::abs(yi) > tolerance; i++) {
        xi -= yi / locally_min_gradient;
        yi = f(xi);
    }
    return xi;
}

} // namespace ccd