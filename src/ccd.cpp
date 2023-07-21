#include "ccd.hpp"
#include "autogen.hpp"

#include <Eigen/Cholesky>

namespace ccd {

namespace {
    static constexpr double EPSILON =
        10 * std::numeric_limits<double>::epsilon();

    std::array<double, 2>
    solve_quadratic_equation(const double a, const double b, const double c)
    {
        if (std::abs(a) < EPSILON) {
            // x = -c / b
            return { { -c / b, -c / b } };
        }

        // x = (-b +- √(b² - 4ac)) / (2a)
        const double discriminant = b * b - 4 * a * c;
        if (discriminant < 0) {
            assert(false);
            return {};
        }
        const double sqrt_discriminant = std::sqrt(discriminant);
        std::array<double, 2> r = { { (-b + sqrt_discriminant) / (2 * a),
                                      (-b - sqrt_discriminant) / (2 * a) } };
        if (r[0] > r[1]) {
            std::swap(r[0], r[1]);
        }
        return r;
    }

    bool is_point_inside_triangle(
        const Eigen::Ref<const Eigen::Vector3d>& p,
        const Eigen::Ref<const Eigen::Vector3d>& t0,
        const Eigen::Ref<const Eigen::Vector3d>& t1,
        const Eigen::Ref<const Eigen::Vector3d>& t2)
    {
        Eigen::Matrix<double, 2, 3> basis;
        basis.row(0) = t1 - t0; // edge 0
        basis.row(1) = t2 - t0; // edge 1
        const Eigen::Matrix2d A = basis * basis.transpose();
        const Eigen::Vector2d b = basis * (p - t0);
        const Eigen::Vector2d x = A.ldlt().solve(b);
        assert((A * x - b).norm() < 1e-10);
        return x[0] >= 0 && x[1] >= 0 && x[0] + x[1] <= 1;
    }

    bool are_edges_intersecting(
        const Eigen::Ref<const Eigen::Vector3d>& ea0,
        const Eigen::Ref<const Eigen::Vector3d>& ea1,
        const Eigen::Ref<const Eigen::Vector3d>& eb0,
        const Eigen::Ref<const Eigen::Vector3d>& eb1)
    {
        const Eigen::Vector3d eb_to_ea = ea0 - eb0;
        const Eigen::Vector3d ea = ea1 - ea0;
        const Eigen::Vector3d eb = eb1 - eb0;

        Eigen::Matrix<double, 2, 2> coefMtr;
        coefMtr(0, 0) = ea.squaredNorm();
        coefMtr(0, 1) = coefMtr(1, 0) = -eb.dot(ea);
        coefMtr(1, 1) = eb.squaredNorm();

        Eigen::Vector2d rhs;
        rhs[0] = -eb_to_ea.dot(ea);
        rhs[1] = eb_to_ea.dot(eb);

        const Eigen::Vector2d x = coefMtr.ldlt().solve(rhs);
        assert((coefMtr * x - rhs).norm() < 1e-10);
        return 0 <= x[0] && x[0] <= 1 && 0 <= x[1] && x[1] <= 1;
    }
} // namespace

std::array<double, 2> CubicEquation::extrema() const
{
    return solve_quadratic_equation(3 * a, 2 * b, c);
}

bool CubicEquation::is_nearly_quadratic() const
{
    return std::abs(a / (b != 0 ? b : 1)) < EPSILON;
}

bool CubicEquation::is_nearly_linear() const
{
    return is_nearly_quadratic() && std::abs(b / (c != 0 ? c : 1)) < EPSILON;
}

bool CubicEquation::is_nearly_constant() const
{
    return is_nearly_linear() && std::abs(c / (d != 0 ? d : 1)) < EPSILON;
}

bool fast_approximate_root_ccd(const CubicEquation d, double& toi)
{
    if (d.is_nearly_constant()) {
        toi = 0;
        return d.d == 0;
    } else if (d.is_nearly_linear()) {
        toi = -d.d / d.c;
        return toi >= 0 && toi <= 1;
    } else if (d.is_nearly_quadratic()) {
        const std::array<double, 2> roots =
            solve_quadratic_equation(d.b, d.c, d.d);
        assert(roots[0] <= roots[1]);
        toi = (0 <= roots[0] && roots[0] <= 1) ? roots[0] : roots[1];
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
    const Eigen::Ref<const Eigen::Vector3d>& p_t0,
    const Eigen::Ref<const Eigen::Vector3d>& t0_t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1_t0,
    const Eigen::Ref<const Eigen::Vector3d>& t2_t0,
    const Eigen::Ref<const Eigen::Vector3d>& p_t1,
    const Eigen::Ref<const Eigen::Vector3d>& t0_t1,
    const Eigen::Ref<const Eigen::Vector3d>& t1_t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2_t1,
    double& toi)
{
    return fast_approximate_root_ccd(
               autogen::point_triangle_ccd_equation(
                   p_t0.x(), p_t0.y(), p_t0.z(), t0_t0.x(), t0_t0.y(),
                   t0_t0.z(), t1_t0.x(), t1_t0.y(), t1_t0.z(), t2_t0.x(),
                   t2_t0.y(), t2_t0.z(), p_t1.x(), p_t1.y(), p_t1.z(),
                   t0_t1.x(), t0_t1.y(), t0_t1.z(), t1_t1.x(), t1_t1.y(),
                   t1_t1.z(), t2_t1.x(), t2_t1.y(), t2_t1.z()),
               toi)
        && is_point_inside_triangle(
               (p_t1 - p_t0) * toi + p_t0, (t0_t1 - t0_t0) * toi + t0_t0,
               (t1_t1 - t1_t0) * toi + t1_t0, (t2_t1 - t2_t0) * toi + t2_t0);
}

bool edge_edge_ccd(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_t0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_t0,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_t0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_t0,
    const Eigen::Ref<const Eigen::Vector3d>& ea0_t1,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_t1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_t1,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_t1,
    double& toi)
{
    return fast_approximate_root_ccd(
               autogen::edge_edge_ccd_equation(
                   ea0_t0.x(), ea0_t0.y(), ea0_t0.z(), ea1_t0.x(), ea1_t0.y(),
                   ea1_t0.z(), eb0_t0.x(), eb0_t0.y(), eb0_t0.z(), eb1_t0.x(),
                   eb1_t0.y(), eb1_t0.z(), ea0_t1.x(), ea0_t1.y(), ea0_t1.z(),
                   ea1_t1.x(), ea1_t1.y(), ea1_t1.z(), eb0_t1.x(), eb0_t1.y(),
                   eb0_t1.z(), eb1_t1.x(), eb1_t1.y(), eb1_t1.z()),
               toi)
        && are_edges_intersecting(
               (ea0_t1 - ea0_t0) * toi + ea0_t0,
               (ea1_t1 - ea1_t0) * toi + ea1_t0,
               (eb0_t1 - eb0_t0) * toi + eb0_t0,
               (eb1_t1 - eb1_t0) * toi + eb1_t0);
}

std::optional<RootInterval>
determine_cubic_root_interval(const CubicEquation& d)
{
    constexpr double t0 = 0, t1 = 1;
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
        return RootInterval(t_max, t_min, d.derivative(-d.b / (3 * d.a)));
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