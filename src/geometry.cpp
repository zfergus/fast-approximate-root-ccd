#include "geometry.hpp"

#include <Eigen/Cholesky>

namespace ccd {

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
} // namespace ccd