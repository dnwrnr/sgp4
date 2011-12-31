#include <Observer.h>
#include <SGP4.h>

#include <iostream>

int main()
{
  Observer obs(51.507406923983446, -0.12773752212524414, 0.05);
  Tle tle = Tle("UK-DMC 2               ", 
    "1 35683U 09041C   11356.17214994  .00000411  00000-0  77254-4 0  6826",
    "2 35683  98.0512 251.8127 0001492  79.4611 280.6776 14.69611889128518");
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
