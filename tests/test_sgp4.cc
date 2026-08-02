#include <gtest/gtest.h>
#include <libsgp4/DecayedException.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/SatelliteException.h>
#include <libsgp4/Tle.h>

using namespace libsgp4;

static const std::string TLE_LINE1_00005 =
    "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
static const std::string TLE_LINE2_00005 =
    "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

static const std::string TLE_LINE1_06251 =
    "1 06251U 62025E   06176.82412014  .00008885  00000-0  12808-3 0  3985";
static const std::string TLE_LINE2_06251 =
    "2 06251  58.0579  54.0425 0030035 139.1568 221.1854 15.56387291  6774";

static const std::string TLE_LINE1_28872 =
    "1 28872U 05037B   05333.02012661  .25992681  00000-0  24476-3 0  1534";
static const std::string TLE_LINE2_28872 =
    "2 28872  96.4736 157.9986 0303955 244.0492 110.6523 16.46015938 10708";

static const std::string TLE_LINE1_88888 =
    "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    87";
static const std::string TLE_LINE2_88888 =
    "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  1058";

TEST(SGP4Basic, ConstructFromTle)
{
    Tle tle(TLE_LINE1_00005, TLE_LINE2_00005);
    SGP4 model(tle);
}

TEST(SGP4Basic, FindPositionAtZero)
{
    Tle tle(TLE_LINE1_00005, TLE_LINE2_00005);
    SGP4 model(tle);
    Eci eci = model.FindPosition(0.0);
    Vector pos = eci.Position();
    Vector vel = eci.Velocity();

    EXPECT_NEAR(pos.x, 7022.46529267, 1e-3);
    EXPECT_NEAR(pos.y, -1400.08296756, 1e-3);
    EXPECT_NEAR(pos.z, 0.03995155, 1e-3);
    EXPECT_NEAR(vel.x, 1.893841015, 1e-6);
    EXPECT_NEAR(vel.y, 6.405893759, 1e-6);
    EXPECT_NEAR(vel.z, 4.534807250, 1e-6);
}

TEST(SGP4Basic, FindPositionAtTsince)
{
    Tle tle(TLE_LINE1_00005, TLE_LINE2_00005);
    SGP4 model(tle);
    Eci eci = model.FindPosition(360.0);
    Vector pos = eci.Position();
    Vector vel = eci.Velocity();

    EXPECT_NEAR(pos.x, -7154.03120202, 1e-3);
    EXPECT_NEAR(pos.y, -3783.17682504, 1e-3);
    EXPECT_NEAR(pos.z, -3536.19412294, 1e-3);
    EXPECT_NEAR(vel.x, 4.741887409, 1e-6);
    EXPECT_NEAR(vel.y, -4.151817765, 1e-6);
    EXPECT_NEAR(vel.z, -2.093935425, 1e-6);
}

TEST(SGP4Basic, FindPositionAtMultipleTimes)
{
    Tle tle(TLE_LINE1_00005, TLE_LINE2_00005);
    SGP4 model(tle);

    double times[] = {0.0, 360.0, 720.0, 1080.0, 1440.0};
    double expectedPos[][3] = {
        {7022.46529267, -1400.08296756, 0.03995155},
        {-7154.03120202, -3783.17682504, -3536.19412294},
        {-7134.59340120, 6531.68641334, 3260.27186483},
        {5568.53901181, 4492.06992591, 3863.87641983},
        {-938.55923943, -6268.18748832, -4294.02924751}
    };

    for (int i = 0; i < 5; i++)
    {
        Eci eci = model.FindPosition(times[i]);
        Vector pos = eci.Position();
        EXPECT_NEAR(pos.x, expectedPos[i][0], 1e-3) << "at t=" << times[i];
        EXPECT_NEAR(pos.y, expectedPos[i][1], 1e-3) << "at t=" << times[i];
        EXPECT_NEAR(pos.z, expectedPos[i][2], 1e-3) << "at t=" << times[i];
    }
}

TEST(SGP4Sat6251, FindPositionAtZero)
{
    Tle tle(TLE_LINE1_06251, TLE_LINE2_06251);
    SGP4 model(tle);
    Eci eci = model.FindPosition(0.0);
    Vector pos = eci.Position();
    Vector vel = eci.Velocity();

    EXPECT_NEAR(pos.x, 3988.31022699, 1e-3);
    EXPECT_NEAR(pos.y, 5498.96657235, 1e-3);
    EXPECT_NEAR(pos.z, 0.90055879, 1e-3);
    EXPECT_NEAR(vel.x, -3.290032738, 1e-6);
    EXPECT_NEAR(vel.y, 2.357652820, 1e-6);
    EXPECT_NEAR(vel.z, 6.496623475, 1e-6);
}

TEST(SGP4Sat6251, FindPositionAt120)
{
    Tle tle(TLE_LINE1_06251, TLE_LINE2_06251);
    SGP4 model(tle);
    Eci eci = model.FindPosition(120.0);
    Vector pos = eci.Position();
    Vector vel = eci.Velocity();

    EXPECT_NEAR(pos.x, -3935.69800083, 1e-3);
    EXPECT_NEAR(pos.y, 409.10980837, 1e-3);
    EXPECT_NEAR(pos.z, 5471.33577327, 1e-3);
    EXPECT_NEAR(vel.x, -3.374784183, 1e-6);
    EXPECT_NEAR(vel.y, -6.635211043, 1e-6);
    EXPECT_NEAR(vel.z, -1.942056221, 1e-6);
}

TEST(SGP4Decayed, ThrowsDecayedException)
{
    Tle tle(TLE_LINE1_28872, TLE_LINE2_28872);
    SGP4 model(tle);

    EXPECT_THROW(model.FindPosition(60.0), DecayedException);
}

TEST(SGP4HighEcc, FindPositionAtZero)
{
    Tle tle(TLE_LINE1_88888, TLE_LINE2_88888);
    SGP4 model(tle);
    Eci eci = model.FindPosition(0.0);
    Vector pos = eci.Position();

    EXPECT_NEAR(pos.x, 2328.9698, 1e-2);
    EXPECT_NEAR(pos.y, -5995.2205, 1e-2);
    EXPECT_NEAR(pos.z, 1719.9730, 1e-2);
}

TEST(SGP4HighEcc, FindPositionAt720)
{
    Tle tle(TLE_LINE1_88888, TLE_LINE2_88888);
    SGP4 model(tle);
    Eci eci = model.FindPosition(720.0);
    Vector pos = eci.Position();

    EXPECT_NEAR(pos.x, 2567.5623, 1e-2);
    EXPECT_NEAR(pos.y, -6112.5038, 1e-2);
    EXPECT_NEAR(pos.z, 713.9637, 1e-2);
}

TEST(SGP4HighEcc, FindPositionAt1440)
{
    Tle tle(TLE_LINE1_88888, TLE_LINE2_88888);
    SGP4 model(tle);
    Eci eci = model.FindPosition(1440.0);
    Vector pos = eci.Position();

    EXPECT_NEAR(pos.x, 2742.5540, 1e-2);
    EXPECT_NEAR(pos.y, -6079.6701, 1e-2);
    EXPECT_NEAR(pos.z, -326.3901, 1e-2);
}

TEST(SGP4SetTle, ChangeModel)
{
    Tle tle1(TLE_LINE1_00005, TLE_LINE2_00005);
    Tle tle2(TLE_LINE1_06251, TLE_LINE2_06251);

    SGP4 model(tle1);
    Eci eci1 = model.FindPosition(0.0);
    Vector pos1 = eci1.Position();

    model.SetTle(tle2);
    Eci eci2 = model.FindPosition(0.0);
    Vector pos2 = eci2.Position();

    EXPECT_NEAR(pos1.x, 7022.46529267, 1e-3);
    EXPECT_NEAR(pos2.x, 3988.31022699, 1e-3);
}

TEST(SGP4Eci, DateTimePosition)
{
    Tle tle(TLE_LINE1_00005, TLE_LINE2_00005);
    SGP4 model(tle);
    DateTime dt = tle.Epoch();
    Eci eci = model.FindPosition(dt);
    Vector pos = eci.Position();

    EXPECT_NEAR(pos.x, 7022.46529267, 1e-3);
    EXPECT_NEAR(pos.y, -1400.08296756, 1e-3);
    EXPECT_NEAR(pos.z, 0.03995155, 1e-3);
}

TEST(SGP4Eci, GetDateTime)
{
    Tle tle(TLE_LINE1_00005, TLE_LINE2_00005);
    SGP4 model(tle);
    Eci eci = model.FindPosition(100.0);
    DateTime dt = eci.GetDateTime();
    EXPECT_EQ(dt.Year(), 2000);
    EXPECT_EQ(dt.Month(), 6);
    EXPECT_EQ(dt.Day(), 27);
}

TEST(SGP4PositionMagnitude, Positive)
{
    Tle tle(TLE_LINE1_00005, TLE_LINE2_00005);
    SGP4 model(tle);
    Eci eci = model.FindPosition(0.0);
    Vector pos = eci.Position();

    double mag = pos.Magnitude();
    EXPECT_GT(mag, 0.0);
    EXPECT_LT(mag, 100000.0);
}

TEST(SGP4VelocityMagnitude, Positive)
{
    Tle tle(TLE_LINE1_00005, TLE_LINE2_00005);
    SGP4 model(tle);
    Eci eci = model.FindPosition(0.0);
    Vector vel = eci.Velocity();

    double mag = vel.Magnitude();
    EXPECT_GT(mag, 0.0);
    EXPECT_LT(mag, 20.0);
}
