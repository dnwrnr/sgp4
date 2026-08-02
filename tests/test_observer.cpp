#include <gtest/gtest.h>
#include <libsgp4/CoordTopocentric.h>
#include <libsgp4/Eci.h>
#include <libsgp4/Observer.h>
#include <libsgp4/SatelliteException.h>
#include <libsgp4/Vector.h>
#include <limits>

using namespace libsgp4;

TEST(ObserverLookAngle, Basic)
{
    DateTime dt(2020, 1, 1, 12, 0, 0);
    Observer observer(40.0, -105.0, 1.6);
    Eci eci(dt, 40.1, -105.2, 500.0);

    CoordTopocentric look = observer.GetLookAngle(eci);
    EXPECT_GE(look.azimuth, 0.0);
    EXPECT_LT(look.azimuth, 2.0 * kPI);
    EXPECT_GE(look.elevation, -kPI / 2.0);
    EXPECT_LE(look.elevation, kPI / 2.0);
    EXPECT_GT(look.range, 0.0);
}

TEST(ObserverLookAngle, ZeroRangeThrows)
{
    DateTime dt(2020, 1, 1, 12, 0, 0);
    Observer observer(0.0, 0.0, 0.0);
    Eci eci(dt, 0.0, 0.0, 0.0);

    EXPECT_THROW(observer.GetLookAngle(eci), SatelliteException);
}

TEST(ObserverLookAngle, NanPositionThrows)
{
    DateTime dt(2020, 1, 1, 12, 0, 0);
    Observer observer(0.0, 0.0, 0.0);
    Vector pos(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0);
    Vector vel(0.0, 0.0, 0.0);
    Eci eci(dt, pos, vel);

    EXPECT_THROW(observer.GetLookAngle(eci), SatelliteException);
}
