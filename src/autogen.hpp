#include "ccd.hpp"

namespace ccd::autogen {

/// @brief Compute the coefficients of the cubic equation for the point-triangle CCD.
/// @return Cubic equation of the point-triangle CCD.
CubicEquation point_triangle_ccd_equation(
    const double x00,
    const double y00,
    const double z00,
    const double x10,
    const double y10,
    const double z10,
    const double x20,
    const double y20,
    const double z20,
    const double x30,
    const double y30,
    const double z30,
    const double x01,
    const double y01,
    const double z01,
    const double x11,
    const double y11,
    const double z11,
    const double x21,
    const double y21,
    const double z21,
    const double x31,
    const double y31,
    const double z31);

/// @brief Compute the coefficients of the cubic equation for the edge-edge CCD.
/// @return Cubic equation of the edge-edge CCD.
CubicEquation edge_edge_ccd_equation(
    const double x00,
    const double y00,
    const double z00,
    const double x10,
    const double y10,
    const double z10,
    const double x20,
    const double y20,
    const double z20,
    const double x30,
    const double y30,
    const double z30,
    const double x01,
    const double y01,
    const double z01,
    const double x11,
    const double y11,
    const double z11,
    const double x21,
    const double y21,
    const double z21,
    const double x31,
    const double y31,
    const double z31);

} // namespace ccd::autogen
