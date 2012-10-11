#include <CoordTopographic.h>
#include <CoordGeodetic.h>
#include <Observer.h>
#include <SGP4.h>

#include <iostream>

int main()
{
    Observer obs(51.507406923983446, -0.12773752212524414, 0.05);
    Tle tle = Tle("MASAT 1                 ",
    "1 25544U 98067A   12285.65009259  .00017228  00000-0  30018-3 0  4501",
    "2 25544 051.6477 262.7396 0017757 155.0745 185.1532 15.50683239796101");
    SGP4 sgp4(tle);

    while (true)
    {
        /*
         * current time
         */
        DateTime now = DateTime::Now(true);
        /*
         * calculate satellite position
         */
        Eci eci = sgp4.FindPosition(now);
        /*
         * get look angle for observer to satellite
         */
        CoordTopographic topo = obs.GetLookAngle(eci);
        /*
         * convert satellite position to geodetic coordinates
         */
        CoordGeodetic geo = eci.ToGeodetic();

        std::cout << now << " " << topo << " " << geo << std::endl;
    };

    return 0;
}
