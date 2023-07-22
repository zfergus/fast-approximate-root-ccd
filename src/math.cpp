#include "math.hpp"

namespace ccd {

std::array<double, 2>
solve_quadratic_equation(const double a, const double b, const double c)
{
    assert(b * b - 4 * a * c >= 0);
    const double tmp = b + sgn(b) * std::sqrt(b * b - 4 * a * c);
    std::array<double, 2> roots = { { -2 * c / tmp, -tmp / (2 * a) } };
    if (roots[0] > roots[1])
        std::swap(roots[0], roots[1]);
    return roots;
}

double
newton_raphson(const CubicEquation& f, const double x0, const double tolerance)
{
    double prev_x, x = x0;
    do {
        prev_x = x;
        x -= std::clamp(f(x) / f.derivative(x), -1.0, 1.0);
    } while (std::abs(x - prev_x) > tolerance);
    return x;
}

double modified_newton_raphson(
    const CubicEquation& f,
    const double x0,
    const double locally_min_gradient,
    const double tolerance)
{
    double prev_x, x = x0;
    do {
        prev_x = x;
        x -= std::clamp(f(x) / locally_min_gradient, -1.0, 1.0);
    } while (std::abs(x - prev_x) > tolerance);
    return x;
}

} // namespace ccd
