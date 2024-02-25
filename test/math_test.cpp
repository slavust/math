#include <gtest/gtest.h>

#include "angle.h"

#include "cylindrical.h"
#include "spherical.h"
#include "polar.h"

#include "vector2.h"
#include "vector3.h"
#include "vector4.h" // just checking that it compiles

#include "matrix2x2.h"
#include "matrix3x3.h"
#include "matrix4x4.h"

#include "euler.h"

#include "matrix_transforms.h"
#include "conv.h"


#include <iostream>

using namespace math;

TEST(Angle, to_radians)
{
	EXPECT_TRUE(std::abs(math::to_radians<double>(180) - math::pi<double>()) < math::eps<double>());
}
TEST(Angle, to_degrees)
{
	EXPECT_TRUE(std::abs(math::to_degrees<double>(math::pi<double>()) - 180) < math::eps<double>());
}

TEST(Cylindrical, Canonize) {
    cylindrical<double> c(1.0, 0.5, 0.0);
    c.canonize();
    EXPECT_DOUBLE_EQ(c.rho, 1.0);
    EXPECT_DOUBLE_EQ(c.phi, 0.5);
    EXPECT_DOUBLE_EQ(c.z, 0.0);

    cylindrical<double> d(-1.0, pi<double>() / 2, 0.0);
    d.canonize();
    EXPECT_DOUBLE_EQ(d.rho, 1.0);
    EXPECT_DOUBLE_EQ(d.phi, -pi<double>() / 2);
    EXPECT_DOUBLE_EQ(d.z, 0.0);

    cylindrical<double> e(-1.0, 3 * pi<double>() / 2, 0.0);
    e.canonize();
    EXPECT_DOUBLE_EQ(e.rho, 1.0);
    EXPECT_DOUBLE_EQ(e.phi, pi<double>() / 2);
    EXPECT_DOUBLE_EQ(e.z, 0.0);

    cylindrical<double> f(0.0, 0.0, 0.0);
    e.canonize();
    EXPECT_DOUBLE_EQ(f.rho, 0.0);
    EXPECT_DOUBLE_EQ(f.phi, 0.0);
    EXPECT_DOUBLE_EQ(f.z, 0.0);
}

TEST(Cylindrical, ToVector3)
{
    const auto validate = [](const vector3<double>& actual, const vector3<double>& expected)
    {
        EXPECT_NEAR(actual.x, expected.x, eps<double>());
        EXPECT_NEAR(actual.y, expected.y, eps<double>());
        EXPECT_NEAR(actual.z, expected.z, eps<double>());
    };


    cylindrical<double> c{2.0, pi<double>() / 2, 3.0};
    vector3<double> expected{0.0, 2.0, 3.0};
    validate(to_vector3(c), expected);
}

TEST(Spherical, Canonize) {
    // Test that the method makes the radius positive
    spherical<double> s1(-1, 0, 0);
    s1.canonize();
    EXPECT_GE(s1.r, 0);

    // Test that the method puts theta in range [0, PI]
    spherical<double> s2(1, 2 * pi<double>(), 0);
    s2.canonize();
    EXPECT_DOUBLE_EQ(s2.theta, 0);

    spherical<double> s3(1, -pi<double>(), 0);
    s3.canonize();
    EXPECT_DOUBLE_EQ(s3.theta, pi<double>());

    // Test that the method puts phi in range [-PI, PI)
    spherical<double> s4(1, 0, 2 * pi<double>());
    s4.canonize();
    EXPECT_DOUBLE_EQ(s4.phi, 0);

    spherical<double> s5(1, 0, -pi<double>());
    s5.canonize();
    EXPECT_DOUBLE_EQ(s5.phi, -pi<double>());

    // Test that the method doesn't change the instance if it's already in the correct range
    spherical<double> s6(1, 0, 0);
    s6.canonize();
    EXPECT_DOUBLE_EQ(s6.r, 1);
    EXPECT_DOUBLE_EQ(s6.theta, 0);
    EXPECT_DOUBLE_EQ(s6.phi, 0);
}

TEST(Spherical, ToVector3)
{
    const auto validate = [](const vector3<double>& actual, const vector3<double>& expected)
    {
        EXPECT_NEAR(actual.x, expected.x, eps<double>());
        EXPECT_NEAR(actual.y, expected.y, eps<double>());
        EXPECT_NEAR(actual.z, expected.z, eps<double>());
    };

    spherical<double> c{2.0, pi<double>() / 2, pi<double>() / 2};
    vector3<double> expected{0.0, 2.0, 0.0};
    validate(to_vector3(c), expected);
}

TEST(Polar, Canonize)
{
    // Test 1: rho = 0, phi = 0
    polar<double> p1(0, 0);
    p1.canonize();
    EXPECT_DOUBLE_EQ(p1.rho, 0);
    EXPECT_DOUBLE_EQ(p1.phi, 0);

    // Test 2: rho = 1, phi = 0
    polar<double> p2(1, 0);
    p2.canonize();
    EXPECT_DOUBLE_EQ(p2.rho, 1);
    EXPECT_DOUBLE_EQ(p2.phi, 0);

    // Test 3: rho = 1, phi = PI
    polar<double> p3(1, pi<double>());
    p3.canonize();
    EXPECT_DOUBLE_EQ(p3.rho, 1);
    EXPECT_DOUBLE_EQ(abs(p3.phi), pi<double>());

    // Test 4: rho = 1, phi = 2*PI
    polar<double> p4(1, 2 * pi<double>());
    p4.canonize();
    EXPECT_DOUBLE_EQ(p4.rho, 1);
    EXPECT_DOUBLE_EQ(p4.phi, 0);

    // Test 5: rho = 1, phi = -2*PI
    polar<double> p6(1, -2 * pi<double>());
    p6.canonize();
    EXPECT_DOUBLE_EQ(p6.rho, 1);
    EXPECT_DOUBLE_EQ(p6.phi, 0);

    // Test 6: rho = -1, phi = 0
    polar<double> p7(-1, 0);
    p7.canonize();
    EXPECT_DOUBLE_EQ(p7.rho, 1);
    EXPECT_DOUBLE_EQ(abs(p7.phi), pi<double>());

    // Test 7: rho = -1, phi = PI
    polar<double> p8(-1, pi<double>());
    p8.canonize();
    EXPECT_DOUBLE_EQ(p8.rho, 1);
    EXPECT_DOUBLE_EQ(p8.phi, 0);

    // Test 8: rho = -1, phi = 2*PI
    polar<double> p9(-1, 2 * pi<double>());
    p9.canonize();
    EXPECT_DOUBLE_EQ(p9.rho, 1);
    EXPECT_DOUBLE_EQ(abs(p9.phi), pi<double>());

    // Test 9: rho = -1, phi = -PI
    polar<double> p10(-1, -pi<double>());
    p10.canonize();
    EXPECT_DOUBLE_EQ(p10.rho, 1);
    EXPECT_DOUBLE_EQ(p10.phi, 0);

    // Test 10: rho = -1, phi = -2*PI
    polar<double> p11(-1, -2 * pi<double>());
    p11.canonize();
    EXPECT_DOUBLE_EQ(p11.rho, 1);
    EXPECT_DOUBLE_EQ(abs(p11.phi), pi<double>());
}

TEST(Polar, ToVector2)
{
    polar<double> p{2.0, pi<double>() / 2};
    vector2<double> v = to_vector2(p);
    vector2<double> expected{0, 2};
    EXPECT_NEAR(v.x, expected.x, eps<double>());
    EXPECT_NEAR(v.y, expected.y, eps<double>());
}

TEST(Vector2, Length) {
    vector2<double> v1(1.0, 2.0);
    EXPECT_DOUBLE_EQ(v1.magnitude(), sqrt(1.0*1.0 + 2.0*2.0));
}

TEST(Vector2, Normalized) {
    vector2<double> v1(1.0, 2.0);
    EXPECT_DOUBLE_EQ(v1.normalized().x, 1.0/sqrt(1.0*1.0 + 2.0*2.0));
    EXPECT_DOUBLE_EQ(v1.normalized().y, 2.0/sqrt(1.0*1.0 + 2.0*2.0));
}

TEST(Vector2, DotProduct) {
    vector2<double> v1(1.0, 2.0);
    vector2<double> v2(3.0, 4.0);
    EXPECT_DOUBLE_EQ(v1.dot(v2), 1.0*3.0 + 2.0*4.0);
}

TEST(Vector3Test, Magnitude) {
    vector3<double> v1(1, 2, 3);
    EXPECT_DOUBLE_EQ(v1.magnitude(), sqrt(14));
}

TEST(Vector3Test, Normalized) {
    vector3<double> v1(1, 2, 3);
    EXPECT_DOUBLE_EQ(v1.normalized().x, 1/sqrt(14));
    EXPECT_DOUBLE_EQ(v1.normalized().y, 2/sqrt(14));
    EXPECT_DOUBLE_EQ(v1.normalized().z, 3/sqrt(14));
}

TEST(Vector3Test, DotProduct) {
    vector3<double> v1(1, 2, 3);
    vector3<double> v2(4, 5, 6);
    EXPECT_DOUBLE_EQ(v1.dot(v2), 1*4 + 2*5 + 3*6);
}

TEST(Vector3Test, CrossProduct) {
    vector3<double> v1(1, 2, 3);
    vector3<double> v2(4, 5, 6);
    EXPECT_DOUBLE_EQ(v1.cross(v2).x, 2*6 - 3*5);
    EXPECT_DOUBLE_EQ(v1.cross(v2).y, 3*4 - 1*6);
    EXPECT_DOUBLE_EQ(v1.cross(v2).z, 1*5 - 2*4);
}


TEST(Matrix2x2, Transposed) {
    matrix2x2<double> m({1, 2, 3, 4});
    matrix2x2<double> expected({1, 3, 2, 4});
    EXPECT_EQ(m.transposed(), expected);
}

TEST(Matrix2x2, Determinant) {
    matrix2x2<double> m({1, 2, 3, 4});
    EXPECT_DOUBLE_EQ(m.determinant(), -2);
}

TEST(Matrix2x2, Adjoint) {
    matrix2x2<double> m(1, 2, 3, 4);
    matrix2x2<double> expected(4, -2, -3, 1);
    matrix2x2<double> adj = m.adjoint();

    for(size_t i = 0; i < 2; ++i)
        for(size_t j = 0; j < 2; ++j)
            EXPECT_DOUBLE_EQ(adj[i][j], expected[i][j]);
}

TEST(Matrix2x2, Inverse) {
    matrix2x2<double> m({1, 2, 3, 4});
    matrix2x2<double> expected({-2, 1, 1.5, -0.5});
    EXPECT_EQ(m.inverse(), expected);
}

TEST(Matrix2x2, Orthonormalize) {
    matrix2x2<double> m({1, 2, 3, 4});
    m.orthonormalize();

    const vector2<double> p{m[0][0], m[0][1]};
    const vector2<double> q{m[1][0], m[1][1]};
    EXPECT_DOUBLE_EQ(p.magnitude(), 1);
    EXPECT_DOUBLE_EQ(q.magnitude(), 1);
    EXPECT_LE(abs(p.dot(q)), eps<double>());
}

TEST(Matrix3x3, Transpose) {
    matrix3x3<double> m{1, 2, 3,
                4, 5, 6,
                7, 8, 9};
    matrix3x3<double> expected{1, 4, 7, 
                               2, 5, 8, 
                               3, 6, 9};
    const auto transposed = m.transposed();

    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(transposed[i][j], expected[i][j]);
}

TEST(Matrix3x3, Multiply) {
    matrix3x3<double> m1{1, 2, 3, 
                        4, 5, 6,
                        7, 8, 9};
    matrix3x3<double> expected{14, 32, 50, 
                               32, 77, 122, 
                               50, 122, 194};

    matrix3x3<double> result = m1 * m1.transposed();
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(result[i][j], expected[i][j]);
}

TEST(Matrix3x3, Inverse) {
    matrix3x3<double> m{0.1, 0., 0.1, 0., 0.2, 0.2, 0., 0., 0.3};
    matrix3x3<double> expected{10., 0., -10./3, 0., 5., -10./3, 0., 0., 10./3};

    matrix3x3<double> inv = m.inverse();
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(inv[i][j], expected[i][j]);
}

TEST(Matrix4x4Test, Transposed) {
    matrix4x4<double> m(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    matrix4x4<double> expected(1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16);
    EXPECT_EQ(m.transposed(), expected);
}

TEST(Matrix4x4Test, Determinant) {
    matrix4x4<double> m(0.1, 0.2, 0.3, 0.4,
                        .5, 0, 0, 0, 
                        0.9, 0.10, 0, 0,
                        0, 0.14, 0.15, 1);
    EXPECT_DOUBLE_EQ(m.determinant(), 0.011999999999999999);
}

TEST(Matrix4x4Test, Inverse) {
    matrix4x4<double> m{0.1, 0., 0.1, 0, 0., 0.2, 0.2, 0, 0., 0., 0.3, 0, 0, 0, 0, 1};
    matrix4x4<double> expected{10., 0., -10./3, 0, 0., 5., -10./3, 0, 0., 0., 10./3, 0, 0, 0, 0, 1};

    const auto inv = m.inverse();

    for(size_t i = 0; i < 4; ++i)
        for(size_t j = 0; j < 4; ++j)
            EXPECT_DOUBLE_EQ(inv[i][j], expected[i][j]);
}

TEST(Matrix4x4Test, Multiplication) {
    matrix4x4<double> m{1, 2, 3, 0,
                        4, 5, 6, 0,
                        7, 8, 9, 0,
                        0, 0, 0, 1};
    matrix4x4<double> expected{14, 32, 50, 0, 
                               32, 77, 122, 0,
                               50, 122, 194, 0,
                               0, 0, 0, 1};
    
    const auto mul = m * m.transposed();
    for(size_t i = 0; i < 4; ++i)
        for(size_t j = 0; j < 4; ++j)
            EXPECT_DOUBLE_EQ(mul[i][j], expected[i][j]);
}


//TEST(EulerTest, CanonizeYaw) {
//    euler<double> euler(1.5 * pi<double>(), 0.0, 0.0);
//    euler.canonize();
//    EXPECT_NEAR(euler.yaw, -0.5 * pi<double>(), eps<double>());
//}
//
//TEST(EulerTest, CanonizePitch) {
//    euler<double> euler(0.0, 1.5 * pi<double>(), 0.0);
//    euler.canonize();
//    EXPECT_NEAR(euler.pitch, -0.5 * pi<double>(), eps<double>());
//}
//
//TEST(EulerTest, CanonizeRoll) {
//    euler<double> euler(0.0, 0.0, 1.5 * pi<double>());
//    euler.canonize();
//    EXPECT_NEAR(euler.roll, -0.5 * pi<double>(), eps<double>());
//}
//
//TEST(EulerTest, CanonizeYawPitch) {
//    euler<double> euler(1.5 * pi<double>(), 1.5 * pi<double>(), 0.0);
//    euler.canonize();
//    EXPECT_NEAR(euler.yaw, -0.5 * pi<double>(), eps<double>());
//    EXPECT_NEAR(euler.pitch, -0.5 * pi<double>(), eps<double>());
//}

//TEST(EulerTest, CanonizeYawRoll) {
//    const auto check_canonized = [](const euler<double>& e)
//    {
//        EXPECT_TRUE(abs(e.yaw) <= pi<double>());
//        EXPECT_TRUE(abs(e.roll) <= pi<double>());
//        EXPECT_TRUE(abs(e.pitch) <= pi<double>() / 2);
//    };
//    const auto validate = [](const vector3<double>& actual, const vector3<double>& expected)
//    {
//        EXPECT_DOUBLE_EQ(actual.x, expected.x);
//        EXPECT_DOUBLE_EQ(actual.y, expected.y);
//        EXPECT_DOUBLE_EQ(actual.z, expected.z);
//    };
//
//    const vector3<double> tested_pos{1, 2, 3};
//
//    euler<double> e(1.5 * pi<double>(), 0.0, 1.5 * pi<double>());
//    const vector3<double> expected_pos = to_rotation_matrix3x3(e) * tested_pos;
//
//    e.canonize();
//    check_canonized(e);
//    validate(expected_pos, to_rotation_matrix3x3(e) * tested_pos);
//}
//
//TEST(EulerTest, CanonizePitcholl) {
//    const auto check_canonized = [](const euler<double>& e)
//    {
//        EXPECT_TRUE(abs(e.yaw) <= pi<double>());
//        EXPECT_TRUE(abs(e.roll) <= pi<double>());
//        EXPECT_TRUE(abs(e.pitch) <= pi<double>() / 2);
//    };
//    const auto validate = [](const vector3<double>& actual, const vector3<double>& expected)
//    {
//        EXPECT_DOUBLE_EQ(actual.x, expected.x);
//        EXPECT_DOUBLE_EQ(actual.y, expected.y);
//        EXPECT_DOUBLE_EQ(actual.z, expected.z);
//    };
//
//    const vector3<double> tested_pos{1, 2, 3};
//
//    euler<double> e(0.0, 1.5 * pi<double>(), 1.5 * pi<double>());
//    const vector3<double> expected_pos = to_rotation_matrix3x3(e) * tested_pos;
//
//    e.canonize();
//    check_canonized(e);
//    validate(expected_pos, to_rotation_matrix3x3(e) * tested_pos);
//}

//TEST(EulerTest, CanonizeYawPitchRoll) {
//    euler<double> euler(1.5 * pi<double>(), 1.5 * pi<double>(), 1.5 * pi<double>());
//    euler.canonize();
//    EXPECT_NEAR(euler.yaw, -0.5 * pi<double>(), eps<double>());
//    EXPECT_NEAR(euler.pitch, -0.5 * pi<double>(), eps<double>());
//    EXPECT_NEAR(euler.roll, -0.5 * pi<double>(), eps<double>());
//}

TEST(MatrixVectorTransformation, Rotation2x2)
{
    const vector2<double> pos{1, 1};

    const vector2<double> expected{-1, -1};
    const vector2<double> r = rotate2D(pi<double>()) * pos;
    EXPECT_DOUBLE_EQ(r.x, expected.x);
    EXPECT_DOUBLE_EQ(r.y, expected.y);
}

TEST(MatrixVectorTransformation, Rotation3x3)
{
    const auto validate = [](const vector3<double>& actual, const vector3<double>& expected)
    {
        EXPECT_DOUBLE_EQ(actual.x, expected.x);
        EXPECT_DOUBLE_EQ(actual.y, expected.y);
        EXPECT_DOUBLE_EQ(actual.z, expected.z);
    };

    const vector3<double> pos{1, 1, 1};
    
    const euler<double> r1{0, pi<double>(), 0};
    const auto p1 = to_rotation_matrix3x3(r1) * pos;
    validate(p1, {1, -1, -1});

    const euler<double> r2(pi<double>(), 0, 0);
    const auto p2 = to_rotation_matrix3x3(r2) * pos;
    validate(p2, {-1, 1, -1});

    const euler<double> r3(0, 0, pi<double>());
    const auto p3 = to_rotation_matrix3x3(r3) * pos;
    validate(p3, {-1, -1, 1});
}

TEST(MatrixVectorTransformation, Rotation4x4)
{
    const auto validate = [](const vector3<double>& actual, const vector3<double>& expected)
    {
        EXPECT_DOUBLE_EQ(actual.x, expected.x);
        EXPECT_DOUBLE_EQ(actual.y, expected.y);
        EXPECT_DOUBLE_EQ(actual.z, expected.z);
    };

    const vector4<double> pos{1, 1, 1, 1};
    
    const euler<double> r1{0, pi<double>(), 0};
    const auto p1 = to_rotation_matrix4x4(r1) * pos;
    validate(p1, {1, -1, -1});

    const euler<double> r2(pi<double>(), 0, 0);
    const auto p2 = to_rotation_matrix4x4(r2) * pos;
    validate(p2, {-1, 1, -1});

    const euler<double> r3(0, 0, pi<double>());
    const auto p3 = to_rotation_matrix4x4(r3) * pos;
    validate(p3, {-1, -1, 1});
}

TEST(MatrixVectorTransformation, Translation3x3)
{
    const auto validate = [](const vector2<double>& actual, const vector2<double>& expected)
    {
        EXPECT_DOUBLE_EQ(actual.x, expected.x);
        EXPECT_DOUBLE_EQ(actual.y, expected.y);
    };

    const vector3<double> pos{1, 1, 1};
    const matrix3x3<double> translation = translate2D<double>(2, 3);
    validate(translation * pos, {3, 4});
}

TEST(MatrixVectorTransformation, Translation4x4)
{
    const auto validate = [](const vector3<double>& actual, const vector3<double>& expected)
    {
        EXPECT_DOUBLE_EQ(actual.x, expected.x);
        EXPECT_DOUBLE_EQ(actual.y, expected.y);
        EXPECT_DOUBLE_EQ(actual.z, expected.z);
    };

    const vector4<double> pos{1, 1, 1, 1};
    const matrix4x4<double> translation = translate3D<double>(1, 2, 3);
    validate(translation * pos, {2, 3, 4});
}