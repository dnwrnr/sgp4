#include "Observer.h"
#include "SGP4.h"

#include <iostream>

int main()
{
  Observer obs(51.507406923983446, -0.12773752212524414, 0.05);
  Tle tle = Tle("ISS (ZARYA)             ",
    "1 25544U 98067A   11146.36888985  .00025753  00000-0  16912-3 0  4201",
    "2 25544  51.6504 272.6534 0003891 329.5510  71.2188 15.75539412717473");
  SGP4 sgp4(tle);

  while(1)
  {
    Julian now;
    Eci eci = sgp4.FindPosition(now);
    CoordTopographic topo = obs.GetLookAngle(eci);
    CoordGeodetic geo = eci.ToGeodetic();
    std::cout << now << " ";
    std::cout << topo << " ";
    std::cout << geo << std::endl;
  };

  return 0;
}
