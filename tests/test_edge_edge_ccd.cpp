#include <catch2/catch_all.hpp>

#include <ccd.hpp>

using namespace ccd;

static const double EPSILON = std::numeric_limits<float>::epsilon();

TEST_CASE("Edge-Edge CCD", "[ccd][edge-edge]")
{
    // e0 = (v0, v1)
    Eigen::Vector3d v0(-1, -1, 0);
    Eigen::Vector3d v1(1, -1, 0);
    // e2 = (v2, v3)
    double e1x = GENERATE(
        -1 - EPSILON, -1, -1 + EPSILON, -0.5, 0, 0.5, 1 - EPSILON, 1,
        1 + EPSILON);
    Eigen::Vector3d v2(e1x, 1, -1);
    Eigen::Vector3d v3(e1x, 1, 1);

    // displacements
    double y_displacement =
        GENERATE(-1.0, 0.0, 1 - EPSILON, 1.0, 1 + EPSILON, 2.0);

    Eigen::Vector3d u0, u1;
    bool is_collision_expected;
    SECTION("moving")
    {
        u0 << 0, y_displacement, 0;
        u1 << 0, -y_displacement, 0;
        is_collision_expected = y_displacement >= 1.0 && e1x >= -1 && e1x <= 1;
    }
    SECTION("fixed")
    {
        u0 << 0, 2 * y_displacement, 0;
        u1.setZero();
        is_collision_expected = y_displacement >= 2.0 && e1x >= -1 && e1x <= 1;
    }

    double toi;
    bool is_colliding =
        edge_edge_ccd(v0, v1, v2, v3, v0 + u0, v1 + u0, v2 + u1, v3 + u1, toi);

    CAPTURE(y_displacement, e1x);
    CHECK(is_colliding >= is_collision_expected);
}

TEST_CASE("Dobule root test case", "[ccd][edge-edge][double-root]")
{
    const Eigen::Vector3d a0s(-3.0022200, 0.2362580, 0.0165247);
    const Eigen::Vector3d a1s(-3.2347850, 0.8312380, -0.1151003);
    const Eigen::Vector3d a0e(-2.8995600, 0.0345838, 0.0638580);
    const Eigen::Vector3d a1e(-3.1716930, 0.6104858, -0.0713340);
    const Eigen::Vector3d b0(-3.0319900, 0.3148750, 0.0000000);
    const Eigen::Vector3d b1(-2.8548800, 0.0900349, 0.0000000);

    bool is_collision_expected = true;

    double toi;
    bool is_colliding = edge_edge_ccd(a0s, a1s, b0, b1, a0e, a1e, b0, b1, toi);

    CHECK(is_colliding >= is_collision_expected);
}

TEST_CASE("Double root test case 2", "[ccd][edge-edge][double-root]")
{
    const Eigen::Vector3d a0s(0, 0, 1);
    const Eigen::Vector3d a1s(0, 1, 1);
    Eigen::Vector3d a0e(1, 1, 0);
    Eigen::Vector3d a1e(0, 0, 0);
    const Eigen::Vector3d b0(0.1, 0.2, 2);
    const Eigen::Vector3d b1(0.1, 0.2, -1);

    double t = GENERATE(0.5, 0.8, 0.88, 0.9, 1.0);
    a0e = (a0e - a0s) * t + a0s;
    a1e = (a1e - a1s) * t + a1s;

    bool is_collision_expected = true;

    double toi;
    bool is_colliding = edge_edge_ccd(a0s, a1s, b0, b1, a0e, a1e, b0, b1, toi);

    CHECK(is_colliding >= is_collision_expected);
}

/*
TEST_CASE("Slow EE CCD", "[ccd][edge-edge]")
{
    Eigen::Vector3d ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1,
        eb1_t1;
    ea0_t0 << 1, 0.50803125, 2.10835646075301e-18;
    ea1_t0 << -2.38233935445388e-18, 0.50803125, 1;
    eb0_t0 << -4.99999999958867e-07, 0.5, 0;
    eb1_t0 << -4.99999999958867e-07, 0.5, 1;
    ea0_t1 << 1, 0.47124375, 4.11078309465837e-18;
    ea1_t1 << -2.8526707189104e-18, 0.47124375, 1;
    eb0_t1 << -4.99999999958867e-07, 0.5, 0;
    eb1_t1 << -4.99999999958867e-07, 0.5, 1;

    // BENCHMARK("compute toi")
    // {
    //     double toi;
    //     bool is_impacting = edge_edge_ccd(
    //         ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
    //         toi);
    // };

    double toi;
    const bool is_impacting = edge_edge_ccd(
        ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, toi);

    CAPTURE(toi);
    CHECK(is_impacting);
}

TEST_CASE("Slow EE CCD 2", "[ccd][edge-edge][slow][thisone]")
{
    Eigen::Vector3d ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1,
        eb1_t1;
    ea0_t0 << 1.00002232466453, 0.500004786049044, -2.06727783590977e-05;
    ea1_t0 << 1.64687846177844e-05, 0.499996645067319, 1.63939999009028e-05;
    eb0_t0 << 1, 0.5, 0;
    eb1_t0 << 0, 0.5, 0;
    ea0_t1 << 1.00294282700155, 0.498652627047143, 0.003626320742036;
    ea1_t1 << -0.00219276550735626, 0.500871179186644, -0.00315828804921928;
    eb0_t1 << 1, 0.5, 0;
    eb1_t1 << 0, 0.5, 0;

    double toi;
    bool is_impacting = edge_edge_ccd(
        ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, toi);

    CAPTURE(toi);
    CHECK(is_impacting);
}
*/