#include "autogen.hpp"

#include <cmath>

namespace ccd::autogen {

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
    const double z31)
{
    const auto t0 = x01 * y11;
    const auto t1 = x01 * y21;
    const auto t2 = x01 * y31;
    const auto t3 = x11 * y01;
    const auto t4 = x11 * y21;
    const auto t5 = x11 * y31;
    const auto t6 = x21 * y01;
    const auto t7 = x21 * y11;
    const auto t8 = x21 * y31;
    const auto t9 = x31 * y01;
    const auto t10 = x31 * y11;
    const auto t11 = x31 * y21;
    const auto t12 = x00 * y10;
    const auto t13 = t12 * z20;
    const auto t14 = x00 * y20;
    const auto t15 = t14 * z30;
    const auto t16 = x00 * y30;
    const auto t17 = t16 * z10;
    const auto t18 = x10 * y00;
    const auto t19 = t18 * z30;
    const auto t20 = x10 * y20;
    const auto t21 = t20 * z00;
    const auto t22 = x10 * y30;
    const auto t23 = t22 * z20;
    const auto t24 = x20 * y00;
    const auto t25 = t24 * z10;
    const auto t26 = x20 * y10;
    const auto t27 = t26 * z30;
    const auto t28 = x20 * y30;
    const auto t29 = t28 * z00;
    const auto t30 = x30 * y00;
    const auto t31 = t30 * z20;
    const auto t32 = x30 * y10;
    const auto t33 = t32 * z00;
    const auto t34 = x30 * y20;
    const auto t35 = t34 * z10;
    const auto t36 = t12 * z30;
    const auto t37 = t14 * z10;
    const auto t38 = t16 * z20;
    const auto t39 = t18 * z20;
    const auto t40 = t20 * z30;
    const auto t41 = t22 * z00;
    const auto t42 = t24 * z30;
    const auto t43 = t26 * z00;
    const auto t44 = t28 * z10;
    const auto t45 = t30 * z10;
    const auto t46 = t32 * z20;
    const auto t47 = t34 * z00;
    const auto t48 = t13 + t15 + t17 + t19 + t21 + t23 + t25 + t27 + t29 + t31
        + t33 + t35 - t36 - t37 - t38 - t39 - t40 - t41 - t42 - t43 - t44 - t45
        - t46 - t47;
    const auto t49 = t12 * z31;
    const auto t50 = x00 * y11;
    const auto t51 = t50 * z30;
    const auto t52 = t14 * z11;
    const auto t53 = x00 * y21;
    const auto t54 = t53 * z10;
    const auto t55 = t16 * z21;
    const auto t56 = x00 * y31;
    const auto t57 = t56 * z20;
    const auto t58 = x01 * y10;
    const auto t59 = t58 * z30;
    const auto t60 = x01 * y20;
    const auto t61 = t60 * z10;
    const auto t62 = x01 * y30;
    const auto t63 = t62 * z20;
    const auto t64 = t18 * z21;
    const auto t65 = x10 * y01;
    const auto t66 = t65 * z20;
    const auto t67 = t20 * z31;
    const auto t68 = x10 * y21;
    const auto t69 = t68 * z30;
    const auto t70 = t22 * z01;
    const auto t71 = x10 * y31;
    const auto t72 = t71 * z00;
    const auto t73 = x11 * y00;
    const auto t74 = t73 * z20;
    const auto t75 = x11 * y20;
    const auto t76 = t75 * z30;
    const auto t77 = x11 * y30;
    const auto t78 = t77 * z00;
    const auto t79 = t24 * z31;
    const auto t80 = x20 * y01;
    const auto t81 = t80 * z30;
    const auto t82 = t26 * z01;
    const auto t83 = x20 * y11;
    const auto t84 = t83 * z00;
    const auto t85 = t28 * z11;
    const auto t86 = x20 * y31;
    const auto t87 = t86 * z10;
    const auto t88 = x21 * y00;
    const auto t89 = t88 * z30;
    const auto t90 = x21 * y10;
    const auto t91 = t90 * z00;
    const auto t92 = x21 * y30;
    const auto t93 = t92 * z10;
    const auto t94 = t30 * z11;
    const auto t95 = x30 * y01;
    const auto t96 = t95 * z10;
    const auto t97 = t32 * z21;
    const auto t98 = x30 * y11;
    const auto t99 = t98 * z20;
    const auto t100 = t34 * z01;
    const auto t101 = x30 * y21;
    const auto t102 = t101 * z00;
    const auto t103 = x31 * y00;
    const auto t104 = t103 * z10;
    const auto t105 = x31 * y10;
    const auto t106 = t105 * z20;
    const auto t107 = x31 * y20;
    const auto t108 = t107 * z00;
    const auto t109 = t100 - t101 * z10 + t102 - t103 * z20 + t104 - t105 * z00
        + t106 - t107 * z10 + t108 - t12 * z21 - t14 * z31 - t16 * z11
        - t18 * z31 - t20 * z01 - t22 * z21 - t24 * z11 - t26 * z31 - t28 * z01
        - t30 * z21 - t32 * z01 - t34 * z11 + t49 - t50 * z20 + t51 + t52
        - t53 * z30 + t54 + t55 - t56 * z10 + t57 - t58 * z20 + t59 - t60 * z30
        + t61 - t62 * z10 + t63 + t64 - t65 * z30 + t66 + t67 - t68 * z00 + t69
        + t70 - t71 * z20 + t72 - t73 * z30 + t74 - t75 * z00 + t76 - t77 * z20
        + t78 + t79 - t80 * z10 + t81 + t82 - t83 * z30 + t84 + t85 - t86 * z00
        + t87 - t88 * z10 + t89 - t90 * z30 + t91 - t92 * z00 + t93 + t94
        - t95 * z20 + t96 + t97 - t98 * z00 + t99;
    const auto t110 = t0 * z20 - t0 * z30 - t1 * z10 + t1 * z30 + t10 * z00
        - t10 * z20 - t101 * z01 + t101 * z11 - t103 * z11 + t103 * z21
        + t105 * z01 - t105 * z21 - t107 * z01 + t107 * z11 - t11 * z00
        + t11 * z10 + t2 * z10 - t2 * z20 - t3 * z20 + t3 * z30 + t4 * z00
        - t4 * z30 - t5 * z00 + t5 * z20 + t50 * z21 - t50 * z31 - t53 * z11
        + t53 * z31 + t56 * z11 - t56 * z21 + t58 * z21 - t58 * z31 + t6 * z10
        - t6 * z30 - t60 * z11 + t60 * z31 + t62 * z11 - t62 * z21 - t65 * z21
        + t65 * z31 + t68 * z01 - t68 * z31 - t7 * z00 + t7 * z30 - t71 * z01
        + t71 * z21 - t73 * z21 + t73 * z31 + t75 * z01 - t75 * z31 - t77 * z01
        + t77 * z21 + t8 * z00 - t8 * z10 + t80 * z11 - t80 * z31 - t83 * z01
        + t83 * z31 + t86 * z01 - t86 * z11 + t88 * z11 - t88 * z31 - t9 * z10
        + t9 * z20 - t90 * z01 + t90 * z31 + t92 * z01 - t92 * z11 - t95 * z11
        + t95 * z21 + t98 * z01 - t98 * z21;
    const auto t111 = 3 * t13 + 3 * t15 + 3 * t17 + 3 * t19 + 3 * t21 + 3 * t23
        + 3 * t25 + 3 * t27 + 3 * t29 + 3 * t31 + 3 * t33 + 3 * t35 - 3 * t36
        - 3 * t37 - 3 * t38 - 3 * t39 - 3 * t40 - 3 * t41 - 3 * t42 - 3 * t43
        - 3 * t44 - 3 * t45 - 3 * t46 - 3 * t47;

    CubicEquation eq;
    eq.a = -t0 * z21 + t0 * z31 + t1 * z11 - t1 * z31 - t10 * z01 + t10 * z21
        + t109 + t11 * z01 - t11 * z11 + t110 - t2 * z11 + t2 * z21 + t3 * z21
        - t3 * z31 - t4 * z01 + t4 * z31 + t48 + t5 * z01 - t5 * z21 - t6 * z11
        + t6 * z31 + t7 * z01 - t7 * z31 - t8 * z01 + t8 * z11 + t9 * z11
        - t9 * z21;
    eq.b = -2 * t100 - 2 * t102 - 2 * t104 - 2 * t106 - 2 * t108 - t110 - t111
        - 2 * t49 - 2 * t51 - 2 * t52 - 2 * t54 - 2 * t55 - 2 * t57 - 2 * t59
        - 2 * t61 - 2 * t63 - 2 * t64 - 2 * t66 - 2 * t67 - 2 * t69 - 2 * t70
        - 2 * t72 - 2 * t74 - 2 * t76 - 2 * t78 - 2 * t79 - 2 * t81 - 2 * t82
        - 2 * t84 - 2 * t85 - 2 * t87 - 2 * t89 - 2 * t91 - 2 * t93 - 2 * t94
        - 2 * t96 - 2 * t97 - 2 * t99 + 2 * x00 * y10 * z21
        + 2 * x00 * y11 * z20 + 2 * x00 * y20 * z31 + 2 * x00 * y21 * z30
        + 2 * x00 * y30 * z11 + 2 * x00 * y31 * z10 + 2 * x01 * y10 * z20
        + 2 * x01 * y20 * z30 + 2 * x01 * y30 * z10 + 2 * x10 * y00 * z31
        + 2 * x10 * y01 * z30 + 2 * x10 * y20 * z01 + 2 * x10 * y21 * z00
        + 2 * x10 * y30 * z21 + 2 * x10 * y31 * z20 + 2 * x11 * y00 * z30
        + 2 * x11 * y20 * z00 + 2 * x11 * y30 * z20 + 2 * x20 * y00 * z11
        + 2 * x20 * y01 * z10 + 2 * x20 * y10 * z31 + 2 * x20 * y11 * z30
        + 2 * x20 * y30 * z01 + 2 * x20 * y31 * z00 + 2 * x21 * y00 * z10
        + 2 * x21 * y10 * z30 + 2 * x21 * y30 * z00 + 2 * x30 * y00 * z21
        + 2 * x30 * y01 * z20 + 2 * x30 * y10 * z01 + 2 * x30 * y11 * z00
        + 2 * x30 * y20 * z11 + 2 * x30 * y21 * z10 + 2 * x31 * y00 * z20
        + 2 * x31 * y10 * z00 + 2 * x31 * y20 * z10;
    eq.c = t109 + t111;
    eq.d = -t48;

    return eq;
}
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
    const double z31)
{
    const auto t0 = x01 * y11;
    const auto t1 = x01 * y21;
    const auto t2 = x01 * y31;
    const auto t3 = x11 * y01;
    const auto t4 = x11 * y21;
    const auto t5 = x11 * y31;
    const auto t6 = x21 * y01;
    const auto t7 = x21 * y11;
    const auto t8 = x21 * y31;
    const auto t9 = x31 * y01;
    const auto t10 = x31 * y11;
    const auto t11 = x31 * y21;
    const auto t12 = x00 * y10;
    const auto t13 = t12 * z20;
    const auto t14 = x00 * y20;
    const auto t15 = t14 * z30;
    const auto t16 = x00 * y30;
    const auto t17 = t16 * z10;
    const auto t18 = x10 * y00;
    const auto t19 = t18 * z30;
    const auto t20 = x10 * y20;
    const auto t21 = t20 * z00;
    const auto t22 = x10 * y30;
    const auto t23 = t22 * z20;
    const auto t24 = x20 * y00;
    const auto t25 = t24 * z10;
    const auto t26 = x20 * y10;
    const auto t27 = t26 * z30;
    const auto t28 = x20 * y30;
    const auto t29 = t28 * z00;
    const auto t30 = x30 * y00;
    const auto t31 = t30 * z20;
    const auto t32 = x30 * y10;
    const auto t33 = t32 * z00;
    const auto t34 = x30 * y20;
    const auto t35 = t34 * z10;
    const auto t36 = t12 * z30;
    const auto t37 = t14 * z10;
    const auto t38 = t16 * z20;
    const auto t39 = t18 * z20;
    const auto t40 = t20 * z30;
    const auto t41 = t22 * z00;
    const auto t42 = t24 * z30;
    const auto t43 = t26 * z00;
    const auto t44 = t28 * z10;
    const auto t45 = t30 * z10;
    const auto t46 = t32 * z20;
    const auto t47 = t34 * z00;
    const auto t48 = t13 + t15 + t17 + t19 + t21 + t23 + t25 + t27 + t29 + t31
        + t33 + t35 - t36 - t37 - t38 - t39 - t40 - t41 - t42 - t43 - t44 - t45
        - t46 - t47;
    const auto t49 = t12 * z31;
    const auto t50 = x00 * y11;
    const auto t51 = t50 * z30;
    const auto t52 = t14 * z11;
    const auto t53 = x00 * y21;
    const auto t54 = t53 * z10;
    const auto t55 = t16 * z21;
    const auto t56 = x00 * y31;
    const auto t57 = t56 * z20;
    const auto t58 = x01 * y10;
    const auto t59 = t58 * z30;
    const auto t60 = x01 * y20;
    const auto t61 = t60 * z10;
    const auto t62 = x01 * y30;
    const auto t63 = t62 * z20;
    const auto t64 = t18 * z21;
    const auto t65 = x10 * y01;
    const auto t66 = t65 * z20;
    const auto t67 = t20 * z31;
    const auto t68 = x10 * y21;
    const auto t69 = t68 * z30;
    const auto t70 = t22 * z01;
    const auto t71 = x10 * y31;
    const auto t72 = t71 * z00;
    const auto t73 = x11 * y00;
    const auto t74 = t73 * z20;
    const auto t75 = x11 * y20;
    const auto t76 = t75 * z30;
    const auto t77 = x11 * y30;
    const auto t78 = t77 * z00;
    const auto t79 = t24 * z31;
    const auto t80 = x20 * y01;
    const auto t81 = t80 * z30;
    const auto t82 = t26 * z01;
    const auto t83 = x20 * y11;
    const auto t84 = t83 * z00;
    const auto t85 = t28 * z11;
    const auto t86 = x20 * y31;
    const auto t87 = t86 * z10;
    const auto t88 = x21 * y00;
    const auto t89 = t88 * z30;
    const auto t90 = x21 * y10;
    const auto t91 = t90 * z00;
    const auto t92 = x21 * y30;
    const auto t93 = t92 * z10;
    const auto t94 = t30 * z11;
    const auto t95 = x30 * y01;
    const auto t96 = t95 * z10;
    const auto t97 = t32 * z21;
    const auto t98 = x30 * y11;
    const auto t99 = t98 * z20;
    const auto t100 = t34 * z01;
    const auto t101 = x30 * y21;
    const auto t102 = t101 * z00;
    const auto t103 = x31 * y00;
    const auto t104 = t103 * z10;
    const auto t105 = x31 * y10;
    const auto t106 = t105 * z20;
    const auto t107 = x31 * y20;
    const auto t108 = t107 * z00;
    const auto t109 = t100 + t102 + t104 + t106 + t108 + t49 + t51 + t52 + t54
        + t55 + t57 + t59 + t61 + t63 + t64 + t66 + t67 + t69 + t70 + t72 + t74
        + t76 + t78 + t79 + t81 + t82 + t84 + t85 + t87 + t89 + t91 + t93 + t94
        + t96 + t97 + t99 - x00 * y10 * z21 - x00 * y11 * z20 - x00 * y20 * z31
        - x00 * y21 * z30 - x00 * y30 * z11 - x00 * y31 * z10 - x01 * y10 * z20
        - x01 * y20 * z30 - x01 * y30 * z10 - x10 * y00 * z31 - x10 * y01 * z30
        - x10 * y20 * z01 - x10 * y21 * z00 - x10 * y30 * z21 - x10 * y31 * z20
        - x11 * y00 * z30 - x11 * y20 * z00 - x11 * y30 * z20 - x20 * y00 * z11
        - x20 * y01 * z10 - x20 * y10 * z31 - x20 * y11 * z30 - x20 * y30 * z01
        - x20 * y31 * z00 - x21 * y00 * z10 - x21 * y10 * z30 - x21 * y30 * z00
        - x30 * y00 * z21 - x30 * y01 * z20 - x30 * y10 * z01 - x30 * y11 * z00
        - x30 * y20 * z11 - x30 * y21 * z10 - x31 * y00 * z20 - x31 * y10 * z00
        - x31 * y20 * z10;
    const auto t110 = t0 * z20 - t0 * z30 - t1 * z10 + t1 * z30 + t10 * z00
        - t10 * z20 - t101 * z01 + t101 * z11 - t103 * z11 + t103 * z21
        + t105 * z01 - t105 * z21 - t107 * z01 + t107 * z11 - t11 * z00
        + t11 * z10 + t2 * z10 - t2 * z20 - t3 * z20 + t3 * z30 + t4 * z00
        - t4 * z30 - t5 * z00 + t5 * z20 + t50 * z21 - t50 * z31 - t53 * z11
        + t53 * z31 + t56 * z11 - t56 * z21 + t58 * z21 - t58 * z31 + t6 * z10
        - t6 * z30 - t60 * z11 + t60 * z31 + t62 * z11 - t62 * z21 - t65 * z21
        + t65 * z31 + t68 * z01 - t68 * z31 - t7 * z00 + t7 * z30 - t71 * z01
        + t71 * z21 - t73 * z21 + t73 * z31 + t75 * z01 - t75 * z31 - t77 * z01
        + t77 * z21 + t8 * z00 - t8 * z10 + t80 * z11 - t80 * z31 - t83 * z01
        + t83 * z31 + t86 * z01 - t86 * z11 + t88 * z11 - t88 * z31 - t9 * z10
        + t9 * z20 - t90 * z01 + t90 * z31 + t92 * z01 - t92 * z11 - t95 * z11
        + t95 * z21 + t98 * z01 - t98 * z21;
    const auto t111 = 3 * t13 + 3 * t15 + 3 * t17 + 3 * t19 + 3 * t21 + 3 * t23
        + 3 * t25 + 3 * t27 + 3 * t29 + 3 * t31 + 3 * t33 + 3 * t35 - 3 * t36
        - 3 * t37 - 3 * t38 - 3 * t39 - 3 * t40 - 3 * t41 - 3 * t42 - 3 * t43
        - 3 * t44 - 3 * t45 - 3 * t46 - 3 * t47;

    CubicEquation eq;
    eq.a = -t0 * z31 - t1 * z11 - t10 * z21 - t109 - t11 * z01 - t110 - t2 * z21
        - t3 * z21 - t4 * z31 - t48 - t5 * z01 - t6 * z31 - t7 * z01 - t8 * z11
        - t9 * z11 + x01 * y11 * z21 + x01 * y21 * z31 + x01 * y31 * z11
        + x11 * y01 * z31 + x11 * y21 * z01 + x11 * y31 * z21 + x21 * y01 * z11
        + x21 * y11 * z31 + x21 * y31 * z01 + x31 * y01 * z21 + x31 * y11 * z01
        + x31 * y21 * z11;
    eq.b = 2 * t100 - 2 * t101 * z10 + 2 * t102 - 2 * t103 * z20 + 2 * t104
        - 2 * t105 * z00 + 2 * t106 - 2 * t107 * z10 + 2 * t108 + t110 + t111
        - 2 * t12 * z21 - 2 * t14 * z31 - 2 * t16 * z11 - 2 * t18 * z31
        - 2 * t20 * z01 - 2 * t22 * z21 - 2 * t24 * z11 - 2 * t26 * z31
        - 2 * t28 * z01 - 2 * t30 * z21 - 2 * t32 * z01 - 2 * t34 * z11
        + 2 * t49 - 2 * t50 * z20 + 2 * t51 + 2 * t52 - 2 * t53 * z30 + 2 * t54
        + 2 * t55 - 2 * t56 * z10 + 2 * t57 - 2 * t58 * z20 + 2 * t59
        - 2 * t60 * z30 + 2 * t61 - 2 * t62 * z10 + 2 * t63 + 2 * t64
        - 2 * t65 * z30 + 2 * t66 + 2 * t67 - 2 * t68 * z00 + 2 * t69 + 2 * t70
        - 2 * t71 * z20 + 2 * t72 - 2 * t73 * z30 + 2 * t74 - 2 * t75 * z00
        + 2 * t76 - 2 * t77 * z20 + 2 * t78 + 2 * t79 - 2 * t80 * z10 + 2 * t81
        + 2 * t82 - 2 * t83 * z30 + 2 * t84 + 2 * t85 - 2 * t86 * z00 + 2 * t87
        - 2 * t88 * z10 + 2 * t89 - 2 * t90 * z30 + 2 * t91 - 2 * t92 * z00
        + 2 * t93 + 2 * t94 - 2 * t95 * z20 + 2 * t96 + 2 * t97 - 2 * t98 * z00
        + 2 * t99;
    eq.c = -t109 - t111;
    eq.d = t48;

    return eq;
}

} // namespace ccd::autogen
