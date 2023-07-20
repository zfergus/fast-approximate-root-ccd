#include "ccd.hpp"
#include "autogen.hpp"

namespace ccd {

static constexpr double EPSILON = 10 * std::numeric_limits<double>::epsilon();

std::array<double, 2>
solve_quadratic_equation(const double a, const double b, const double c)
{
    if (std::abs(a) < EPSILON) {
        // x = -c / b
        return { -c / b, -c / b };
    }

    // x = (-b +- √(b² - 4ac)) / (2a)
    const double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        assert(false);
        return {};
    }
    const double sqrt_discriminant = std::sqrt(discriminant);
    std::array<double, 2> r = { (-b + sqrt_discriminant) / (2 * a),
                                (-b - sqrt_discriminant) / (2 * a) };
    if (r[0] > r[1]) {
        std::swap(r[0], r[1]);
    }
    return r;
}

std::array<double, 2> CubicEquation::extrema() const
{
    return solve_quadratic_equation(3 * a, 2 * b, c);
}

bool fast_approximate_root_ccd(const CubicEquation d, double& toi)
{
    if (d.a == 0) {
        if (d.b == 0) {
            toi = d.c == 0 ? d.d : (-d.d / d.c);
        } else {
            const std::array<double, 2> roots =
                solve_quadratic_equation(d.b, d.c, d.d);
            assert(roots[0] <= roots[1]);
            toi = (0 <= roots[0] && roots[0] <= 1) ? roots[0] : roots[1];
        }
        return toi >= 0 && toi <= 1;
    }

    const auto root_interval = determine_cubic_root_interval(d);
    if (!root_interval || root_interval->a > root_interval->b) {
        return false;
    }

    toi = modified_newton_raphson(
        d, root_interval->a, root_interval->min_gradient);

    return toi >= 0 && toi <= 1;
}

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
    return fast_approximate_root_ccd(
        autogen::point_triangle_ccd_equation(
            p_t0.x(), p_t0.y(), p_t0.z(), t0_t0.x(), t0_t0.y(), t0_t0.z(),
            t1_t0.x(), t1_t0.y(), t1_t0.z(), t2_t0.x(), t2_t0.y(), t2_t0.z(),
            p_t1.x(), p_t1.y(), p_t1.z(), t0_t1.x(), t0_t1.y(), t0_t1.z(),
            t1_t1.x(), t1_t1.y(), t1_t1.z(), t2_t1.x(), t2_t1.y(), t2_t1.z()),
        toi);
}

bool edge_edge_ccd(
    const Eigen::Vector3d& ea0_t0,
    const Eigen::Vector3d& ea1_t0,
    const Eigen::Vector3d& eb0_t0,
    const Eigen::Vector3d& eb1_t0,
    const Eigen::Vector3d& ea0_t1,
    const Eigen::Vector3d& ea1_t1,
    const Eigen::Vector3d& eb0_t1,
    const Eigen::Vector3d& eb1_t1,
    double& toi)
{
    return fast_approximate_root_ccd(
        autogen::edge_edge_ccd_equation(
            ea0_t0.x(), ea0_t0.y(), ea0_t0.z(), ea1_t0.x(), ea1_t0.y(),
            ea1_t0.z(), eb0_t0.x(), eb0_t0.y(), eb0_t0.z(), eb1_t0.x(),
            eb1_t0.y(), eb1_t0.z(), ea0_t1.x(), ea0_t1.y(), ea0_t1.z(),
            ea1_t1.x(), ea1_t1.y(), ea1_t1.z(), eb0_t1.x(), eb0_t1.y(),
            eb0_t1.z(), eb1_t1.x(), eb1_t1.y(), eb1_t1.z()),
        toi);
}

std::optional<RootInterval>
determine_cubic_root_interval(const CubicEquation& d)
{
    double t0 = 0, t1 = 1;
    assert(d(t0) != 0);

    const auto [tm0, tm1] = d.extrema();
    assert(tm0 <= tm1);

    assert(d.a != 0);
    if (d.a > 0) {
        // Case 1 (Interval II)
        const double t_min = tm1, t_max = tm0;
        if (d(t_min) > 0 || t0 >= t_min) {
            // assert(d(t0) > 0);
            return std::nullopt;
        }
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