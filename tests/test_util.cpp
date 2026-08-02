#include <gtest/gtest.h>
#include <libsgp4/Globals.h>
#include <libsgp4/Util.h>
#include <limits>

using namespace libsgp4;

TEST(UtilClamp, WithinRange)
{
    EXPECT_DOUBLE_EQ(Util::Clamp(0.0, -1.0, 1.0), 0.0);
    EXPECT_DOUBLE_EQ(Util::Clamp(0.5, -1.0, 1.0), 0.5);
    EXPECT_DOUBLE_EQ(Util::Clamp(-0.5, -1.0, 1.0), -0.5);
}

TEST(UtilClamp, AtBounds)
{
    EXPECT_DOUBLE_EQ(Util::Clamp(-1.0, -1.0, 1.0), -1.0);
    EXPECT_DOUBLE_EQ(Util::Clamp(1.0, -1.0, 1.0), 1.0);
}

TEST(UtilClamp, BelowMin)
{
    EXPECT_DOUBLE_EQ(Util::Clamp(-2.0, -1.0, 1.0), -1.0);
}

TEST(UtilClamp, AboveMax)
{
    EXPECT_DOUBLE_EQ(Util::Clamp(2.0, -1.0, 1.0), 1.0);
}

TEST(UtilClamp, NaN)
{
    EXPECT_DOUBLE_EQ(Util::Clamp(std::numeric_limits<double>::quiet_NaN(), -1.0, 1.0), -1.0);
}

TEST(UtilMod, PositiveDividend)
{
    EXPECT_DOUBLE_EQ(Util::Mod(7.0, 3.0), 1.0);
    EXPECT_DOUBLE_EQ(Util::Mod(10.0, 5.0), 0.0);
    EXPECT_DOUBLE_EQ(Util::Mod(1.0, 3.0), 1.0);
}

TEST(UtilMod, NegativeDividend)
{
    EXPECT_DOUBLE_EQ(Util::Mod(-3.0, 4.0), 1.0);
    EXPECT_DOUBLE_EQ(Util::Mod(-1.0, 3.0), 2.0);
    EXPECT_DOUBLE_EQ(Util::Mod(-7.0, 3.0), 2.0);
}

TEST(UtilMod, ZeroDivisor)
{
    EXPECT_DOUBLE_EQ(Util::Mod(5.0, 0.0), 5.0);
    EXPECT_DOUBLE_EQ(Util::Mod(-5.0, 0.0), -5.0);
}

TEST(UtilMod, LargeValues)
{
    EXPECT_DOUBLE_EQ(Util::Mod(360.0, 360.0), 0.0);
    EXPECT_DOUBLE_EQ(Util::Mod(720.0, 360.0), 0.0);
    EXPECT_DOUBLE_EQ(Util::Mod(361.0, 360.0), 1.0);
}

TEST(UtilWrapNegPosPI, Zero)
{
    EXPECT_NEAR(Util::WrapNegPosPI(0.0), 0.0, 1e-10);
}

TEST(UtilWrapNegPosPI, WithinRange)
{
    EXPECT_NEAR(Util::WrapNegPosPI(1.0), 1.0, 1e-10);
    EXPECT_NEAR(Util::WrapNegPosPI(-1.0), -1.0, 1e-10);
}

TEST(UtilWrapNegPosPI, AtPI)
{
    EXPECT_NEAR(Util::WrapNegPosPI(kPI), -kPI, 1e-10);
    EXPECT_NEAR(Util::WrapNegPosPI(-kPI), -kPI, 1e-10);
}

TEST(UtilWrapNegPosPI, BeyondPI)
{
    EXPECT_NEAR(Util::WrapNegPosPI(2.0 * kPI), 0.0, 1e-10);
    EXPECT_NEAR(Util::WrapNegPosPI(-2.0 * kPI), 0.0, 1e-10);
    EXPECT_NEAR(Util::WrapNegPosPI(3.0 * kPI), -kPI, 1e-10);
    EXPECT_NEAR(Util::WrapNegPosPI(-3.0 * kPI), -kPI, 1e-10);
}

TEST(UtilWrapTwoPI, Zero)
{
    EXPECT_NEAR(Util::WrapTwoPI(0.0), 0.0, 1e-10);
}

TEST(UtilWrapTwoPI, WithinRange)
{
    EXPECT_NEAR(Util::WrapTwoPI(1.0), 1.0, 1e-10);
}

TEST(UtilWrapTwoPI, AtTwoPI)
{
    EXPECT_NEAR(Util::WrapTwoPI(kTWOPI), 0.0, 1e-10);
}

TEST(UtilWrapTwoPI, BeyondTwoPI)
{
    EXPECT_NEAR(Util::WrapTwoPI(2.0 * kTWOPI), 0.0, 1e-10);
    EXPECT_NEAR(Util::WrapTwoPI(3.0 * kTWOPI), 0.0, 1e-10);
}

TEST(UtilWrapTwoPI, Negative)
{
    EXPECT_NEAR(Util::WrapTwoPI(-1.0), kTWOPI - 1.0, 1e-10);
    EXPECT_NEAR(Util::WrapTwoPI(-kTWOPI), 0.0, 1e-10);
}

TEST(UtilWrapNegPos180, Zero)
{
    EXPECT_DOUBLE_EQ(Util::WrapNegPos180(0.0), 0.0);
}

TEST(UtilWrapNegPos180, WithinRange)
{
    EXPECT_DOUBLE_EQ(Util::WrapNegPos180(90.0), 90.0);
    EXPECT_DOUBLE_EQ(Util::WrapNegPos180(-90.0), -90.0);
}

TEST(UtilWrapNegPos180, Beyond180)
{
    EXPECT_DOUBLE_EQ(Util::WrapNegPos180(360.0), 0.0);
    EXPECT_DOUBLE_EQ(Util::WrapNegPos180(-360.0), 0.0);
    EXPECT_DOUBLE_EQ(Util::WrapNegPos180(180.0), -180.0);
}

TEST(UtilWrap360, Zero)
{
    EXPECT_DOUBLE_EQ(Util::Wrap360(0.0), 0.0);
}

TEST(UtilWrap360, WithinRange)
{
    EXPECT_DOUBLE_EQ(Util::Wrap360(180.0), 180.0);
}

TEST(UtilWrap360, Beyond360)
{
    EXPECT_DOUBLE_EQ(Util::Wrap360(360.0), 0.0);
    EXPECT_DOUBLE_EQ(Util::Wrap360(720.0), 0.0);
    EXPECT_DOUBLE_EQ(Util::Wrap360(361.0), 1.0);
}

TEST(UtilWrap360, Negative)
{
    EXPECT_DOUBLE_EQ(Util::Wrap360(-1.0), 359.0);
    EXPECT_DOUBLE_EQ(Util::Wrap360(-360.0), 0.0);
}

TEST(UtilDegreesToRadians, Zero)
{
    EXPECT_DOUBLE_EQ(Util::DegreesToRadians(0.0), 0.0);
}

TEST(UtilDegreesToRadians, KnownValues)
{
    EXPECT_NEAR(Util::DegreesToRadians(90.0), kPI / 2.0, 1e-10);
    EXPECT_NEAR(Util::DegreesToRadians(180.0), kPI, 1e-10);
    EXPECT_NEAR(Util::DegreesToRadians(270.0), 3.0 * kPI / 2.0, 1e-10);
    EXPECT_NEAR(Util::DegreesToRadians(360.0), kTWOPI, 1e-10);
}

TEST(UtilDegreesToRadians, Negative)
{
    EXPECT_NEAR(Util::DegreesToRadians(-90.0), -kPI / 2.0, 1e-10);
    EXPECT_NEAR(Util::DegreesToRadians(-180.0), -kPI, 1e-10);
}

TEST(UtilRadiansToDegrees, Zero)
{
    EXPECT_DOUBLE_EQ(Util::RadiansToDegrees(0.0), 0.0);
}

TEST(UtilRadiansToDegrees, KnownValues)
{
    EXPECT_NEAR(Util::RadiansToDegrees(kPI / 2.0), 90.0, 1e-10);
    EXPECT_NEAR(Util::RadiansToDegrees(kPI), 180.0, 1e-10);
    EXPECT_NEAR(Util::RadiansToDegrees(3.0 * kPI / 2.0), 270.0, 1e-10);
    EXPECT_NEAR(Util::RadiansToDegrees(kTWOPI), 360.0, 1e-10);
}

TEST(UtilDegreesRadiansRoundTrip, Values)
{
    EXPECT_NEAR(Util::RadiansToDegrees(Util::DegreesToRadians(45.0)), 45.0, 1e-10);
    EXPECT_NEAR(Util::DegreesToRadians(Util::RadiansToDegrees(kPI / 3.0)), kPI / 3.0, 1e-10);
}

TEST(UtilAcTan, Quadrant1)
{
    EXPECT_NEAR(Util::AcTan(1.0, 1.0), kPI / 4.0, 1e-10);
    EXPECT_NEAR(Util::AcTan(0.0, 1.0), 0.0, 1e-10);
}

TEST(UtilAcTan, Quadrant2)
{
    EXPECT_NEAR(Util::AcTan(1.0, -1.0), 3.0 * kPI / 4.0, 1e-10);
}

TEST(UtilAcTan, Quadrant3)
{
    EXPECT_NEAR(Util::AcTan(-1.0, -1.0), 5.0 * kPI / 4.0, 1e-10);
}

TEST(UtilAcTan, Quadrant4)
{
    EXPECT_NEAR(Util::AcTan(-1.0, 1.0), -kPI / 4.0, 1e-10);
}

TEST(UtilAcTan, CosxZero)
{
    EXPECT_NEAR(Util::AcTan(1.0, 0.0), kPI / 2.0, 1e-10);
    EXPECT_NEAR(Util::AcTan(-1.0, 0.0), 3.0 * kPI / 2.0, 1e-10);
}

TEST(UtilTrim, LeadingSpaces)
{
    std::string s = "  hello";
    Util::TrimLeft(s);
    EXPECT_EQ(s, "hello");
}

TEST(UtilTrim, TrailingSpaces)
{
    std::string s = "hello  ";
    Util::TrimRight(s);
    EXPECT_EQ(s, "hello");
}

TEST(UtilTrim, BothSpaces)
{
    std::string s = "  hello  ";
    Util::Trim(s);
    EXPECT_EQ(s, "hello");
}

TEST(UtilTrim, NoSpaces)
{
    std::string s = "hello";
    Util::Trim(s);
    EXPECT_EQ(s, "hello");
}

TEST(UtilTrim, AllSpaces)
{
    std::string s = "   ";
    Util::Trim(s);
    EXPECT_EQ(s, "");
}

TEST(UtilTrim, TabsAndNewlines)
{
    std::string s = "\t\nhello\t\n";
    Util::Trim(s);
    EXPECT_EQ(s, "hello");
}

TEST(UtilFromString, ValidInt)
{
    int val;
    EXPECT_TRUE(Util::FromString("42", val));
    EXPECT_EQ(val, 42);
}

TEST(UtilFromString, ValidDouble)
{
    double val;
    EXPECT_TRUE(Util::FromString("3.14", val));
    EXPECT_DOUBLE_EQ(val, 3.14);
}

TEST(UtilFromString, InvalidString)
{
    double val;
    EXPECT_FALSE(Util::FromString("abc", val));
}
