#include <CoordTopographic.h>
#include <CoordGeodetic.h>
#include <Observer.h>
#include <SGP4.h>

#include <iostream>

int main()
{
    Observer obs(51.507406923983446, -0.12773752212524414, 0.05);
    Tle tle = Tle("UK-DMC 2                ",
        "1 35683U 09041C   12289.23158813  .00000484  00000-0  89219-4 0  5863",
        "2 35683  98.0221 185.3682 0001499 100.5295 259.6088 14.69819587172294");
    SGP4 sgp4(tle);

    std::cout << tle << std::endl;

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
