#include <gtest/gtest.h>
#include <libsgp4/Vector.h>

using namespace libsgp4;

TEST(VectorConstruction, Default)
{
    Vector v;
    EXPECT_DOUBLE_EQ(v.x, 0.0);
    EXPECT_DOUBLE_EQ(v.y, 0.0);
    EXPECT_DOUBLE_EQ(v.z, 0.0);
    EXPECT_DOUBLE_EQ(v.w, 0.0);
}

TEST(VectorConstruction, ThreeArgs)
{
    Vector v(1.0, 2.0, 3.0);
    EXPECT_DOUBLE_EQ(v.x, 1.0);
    EXPECT_DOUBLE_EQ(v.y, 2.0);
    EXPECT_DOUBLE_EQ(v.z, 3.0);
    EXPECT_DOUBLE_EQ(v.w, 0.0);
}

TEST(VectorConstruction, FourArgs)
{
    Vector v(1.0, 2.0, 3.0, 4.0);
    EXPECT_DOUBLE_EQ(v.x, 1.0);
    EXPECT_DOUBLE_EQ(v.y, 2.0);
    EXPECT_DOUBLE_EQ(v.z, 3.0);
    EXPECT_DOUBLE_EQ(v.w, 4.0);
}

TEST(VectorConstruction, Copy)
{
    Vector v1(1.0, 2.0, 3.0, 4.0);
    Vector v2(v1);
    EXPECT_DOUBLE_EQ(v2.x, 1.0);
    EXPECT_DOUBLE_EQ(v2.y, 2.0);
    EXPECT_DOUBLE_EQ(v2.z, 3.0);
    EXPECT_DOUBLE_EQ(v2.w, 4.0);
}

TEST(VectorAssignment, Equal)
{
    Vector v1(1.0, 2.0, 3.0, 4.0);
    Vector v2;
    v2 = v1;
    EXPECT_DOUBLE_EQ(v2.x, 1.0);
    EXPECT_DOUBLE_EQ(v2.y, 2.0);
    EXPECT_DOUBLE_EQ(v2.z, 3.0);
    EXPECT_DOUBLE_EQ(v2.w, 4.0);
}

TEST(VectorMagnitude, ZeroVector)
{
    Vector v;
    EXPECT_DOUBLE_EQ(v.Magnitude(), 0.0);
}

TEST(VectorMagnitude, UnitVector)
{
    Vector v(1.0, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(v.Magnitude(), 1.0);

    Vector v2(0.0, 1.0, 0.0);
    EXPECT_DOUBLE_EQ(v2.Magnitude(), 1.0);

    Vector v3(0.0, 0.0, 1.0);
    EXPECT_DOUBLE_EQ(v3.Magnitude(), 1.0);
}

TEST(VectorMagnitude, KnownValues)
{
    Vector v(3.0, 4.0, 0.0);
    EXPECT_DOUBLE_EQ(v.Magnitude(), 5.0);

    Vector v2(1.0, 2.0, 2.0);
    EXPECT_DOUBLE_EQ(v2.Magnitude(), 3.0);
}

TEST(VectorMagnitude, NegativeComponents)
{
    Vector v(-3.0, -4.0, 0.0);
    EXPECT_DOUBLE_EQ(v.Magnitude(), 5.0);
}

TEST(VectorMagnitude, IgnoresW)
{
    Vector v(1.0, 0.0, 0.0, 999.0);
    EXPECT_DOUBLE_EQ(v.Magnitude(), 1.0);
}

TEST(VectorDot, Orthogonal)
{
    Vector v1(1.0, 0.0, 0.0);
    Vector v2(0.0, 1.0, 0.0);
    EXPECT_DOUBLE_EQ(v1.Dot(v2), 0.0);
}

TEST(VectorDot, Parallel)
{
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(1.0, 2.0, 3.0);
    EXPECT_DOUBLE_EQ(v1.Dot(v2), 14.0);
}

TEST(VectorDot, KnownValues)
{
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(4.0, 5.0, 6.0);
    EXPECT_DOUBLE_EQ(v1.Dot(v2), 32.0);
}

TEST(VectorDot, IgnoresW)
{
    Vector v1(1.0, 0.0, 0.0, 999.0);
    Vector v2(1.0, 0.0, 0.0, 999.0);
    EXPECT_DOUBLE_EQ(v1.Dot(v2), 1.0);
}

TEST(VectorSubtract, Basic)
{
    Vector v1(5.0, 5.0, 5.0);
    Vector v2(1.0, 2.0, 3.0);
    Vector result = v1 - v2;
    EXPECT_DOUBLE_EQ(result.x, 4.0);
    EXPECT_DOUBLE_EQ(result.y, 3.0);
    EXPECT_DOUBLE_EQ(result.z, 2.0);
}

TEST(VectorSubtract, WIsZero)
{
    Vector v1(1.0, 1.0, 1.0, 5.0);
    Vector v2(0.0, 0.0, 0.0, 3.0);
    Vector result = v1 - v2;
    EXPECT_DOUBLE_EQ(result.w, 0.0);
}

TEST(VectorSubtract, IdenticalVectors)
{
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(1.0, 2.0, 3.0);
    Vector result = v1 - v2;
    EXPECT_DOUBLE_EQ(result.x, 0.0);
    EXPECT_DOUBLE_EQ(result.y, 0.0);
    EXPECT_DOUBLE_EQ(result.z, 0.0);
}

TEST(VectorToString, Format)
{
    Vector v(1.0, 2.0, 3.0, 4.0);
    std::string str = v.ToString();
    EXPECT_NE(str.find("1.000"), std::string::npos);
    EXPECT_NE(str.find("2.000"), std::string::npos);
    EXPECT_NE(str.find("3.000"), std::string::npos);
    EXPECT_NE(str.find("4.000"), std::string::npos);
}
