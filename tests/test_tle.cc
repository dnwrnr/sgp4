#include <gtest/gtest.h>
#include <libsgp4/Tle.h>
#include <libsgp4/TleException.h>

using namespace libsgp4;

static const std::string VALID_LINE1 =
    "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
static const std::string VALID_LINE2 =
    "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

TEST(TleConstruction, ValidTwoLine)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_EQ(tle.NoradNumber(), 5u);
}

TEST(TleConstruction, ValidWithName)
{
    Tle tle("TEST SAT", VALID_LINE1, VALID_LINE2);
    EXPECT_EQ(tle.Name(), "TEST SAT");
    EXPECT_EQ(tle.NoradNumber(), 5u);
}

TEST(TleConstruction, InvalidLineLengthLine1)
{
    std::string shortLine = "1 00005";
    EXPECT_THROW(Tle(shortLine, VALID_LINE2), TleException);
}

TEST(TleConstruction, InvalidLineLengthLine2)
{
    std::string shortLine = "2 00005";
    EXPECT_THROW(Tle(VALID_LINE1, shortLine), TleException);
}

TEST(TleConstruction, WrongLine1Prefix)
{
    std::string badLine = VALID_LINE1;
    badLine[0] = '2';
    EXPECT_THROW(Tle(badLine, VALID_LINE2), TleException);
}

TEST(TleConstruction, WrongLine2Prefix)
{
    std::string badLine = VALID_LINE2;
    badLine[0] = '1';
    EXPECT_THROW(Tle(VALID_LINE1, badLine), TleException);
}

TEST(TleConstruction, NORADNumberMismatch)
{
    std::string line2Different = VALID_LINE2;
    line2Different[2] = '9';
    EXPECT_THROW(Tle(VALID_LINE1, line2Different), TleException);
}

TEST(TleEpoch, YearBefore57)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    DateTime epoch = tle.Epoch();
    EXPECT_EQ(epoch.Year(), 2000);
}

TEST(TleEpoch, YearAfter57)
{
    std::string line1 = "1 04632U 70093B   04031.91070959 -.00000084  00000-0  10000-3 0  9955";
    std::string line2 = "2 04632  11.4628 273.1101 1450506 207.6000 143.9350  1.20231981 44145";
    Tle tle(line1, line2);
    DateTime epoch = tle.Epoch();
    EXPECT_EQ(epoch.Year(), 2004);
}

TEST(TleEpoch, DayOfYear)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    DateTime epoch = tle.Epoch();
    EXPECT_EQ(epoch.Month(), 6);
    EXPECT_EQ(epoch.Day(), 27);
}

TEST(TleFields, Inclination)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_NEAR(tle.Inclination(true), 34.2682, 1e-4);
}

TEST(TleFields, RightAscendingNode)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_NEAR(tle.RightAscendingNode(true), 348.7242, 1e-4);
}

TEST(TleFields, Eccentricity)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_NEAR(tle.Eccentricity(), 0.1859667, 1e-7);
}

TEST(TleFields, ArgumentPerigee)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_NEAR(tle.ArgumentPerigee(true), 331.7664, 1e-4);
}

TEST(TleFields, MeanAnomaly)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_NEAR(tle.MeanAnomaly(true), 19.3264, 1e-4);
}

TEST(TleFields, MeanMotion)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_NEAR(tle.MeanMotion(), 10.82419157, 1e-4);
}

TEST(TleFields, OrbitNumber)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_EQ(tle.OrbitNumber(), 41366u);
}

TEST(TleFields, MeanMotionDt2)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_NEAR(tle.MeanMotionDt2(), 2.3e-7, 1e-10);
}

TEST(TleFields, MeanMotionDdt6)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_NEAR(tle.MeanMotionDdt6(), 0.0, 1e-10);
}

TEST(TleFields, BStar)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_NEAR(tle.BStar(), 2.8098e-5, 1e-8);
}

TEST(TleFields, IntDesignator)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    EXPECT_EQ(tle.IntDesignator(), "58002B  ");
}

TEST(TleLineLength, Expected)
{
    EXPECT_EQ(Tle::LineLength(), 69u);
    EXPECT_EQ(VALID_LINE1.length(), 69u);
    EXPECT_EQ(VALID_LINE2.length(), 69u);
}

TEST(TleMultipleTLEs, IndependentParsing)
{
    std::string line1A = "1 04632U 70093B   04031.91070959 -.00000084  00000-0  10000-3 0  9955";
    std::string line2A = "2 04632  11.4628 273.1101 1450506 207.6000 143.9350  1.20231981 44145";

    Tle tle1(VALID_LINE1, VALID_LINE2);
    Tle tle2(line1A, line2A);

    EXPECT_EQ(tle1.NoradNumber(), 5u);
    EXPECT_EQ(tle2.NoradNumber(), 4632u);
    EXPECT_NEAR(tle1.Inclination(true), 34.2682, 1e-4);
    EXPECT_NEAR(tle2.Inclination(true), 11.4628, 1e-4);
}

TEST(TleNegativeBSTAR, NegativeValue)
{
    std::string line1 = "1 21897U 92011A   06176.02341244 -.00001273  00000-0 -13525-3 0  3044";
    std::string line2 = "2 21897  62.1749 198.0096 7421690 253.0462  20.1561  2.01269994104880";
    Tle tle(line1, line2);
    EXPECT_LT(tle.BStar(), 0.0);
}

TEST(TleMeanMotionDdt6, PositiveExponent)
{
    std::string line1 = "1 11801U          80230.29629788  .01431103  00000-0  14311-1      13";
    std::string line2 = "2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848    13";
    Tle tle(line1, line2);
    EXPECT_NEAR(tle.MeanMotionDdt6(), 0.0, 1e-10);
    EXPECT_NEAR(tle.BStar(), 0.014311, 1e-6);
}

TEST(TleToString, Output)
{
    Tle tle(VALID_LINE1, VALID_LINE2);
    std::string str = tle.ToString();
    EXPECT_NE(str.find("Norad Number"), std::string::npos);
}

TEST(TleMeanMotion, HighEcc)
{
    std::string line1 = "1 23333U 94071A   94305.49999999 -.00172956  26967-3  10000-3 0    15";
    std::string line2 = "2 23333  28.7490   2.3720 9728298  30.4360   1.3500  0.07309491    70";
    Tle tle(line1, line2);
    EXPECT_NEAR(tle.MeanMotion(), 0.07309491, 1e-4);
    EXPECT_NEAR(tle.Eccentricity(), 0.9728298, 1e-7);
}
