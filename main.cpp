#include "Julian.h"
#include "Tle.h"
#include "SGDP4.h"
#include "Globals.h"
#include "Observer.h"
#include "Coord.h"

#include <list>
#include <string>
#include <iomanip>

void FindSatellite(const Julian& time_start, const Julian& time_end) {

    /*
     * half a second
     */
    static const double delta = 1.0 / (60.0 * 60.0 * 24.0 * 2.0);

    //while (fabs(time_end - time_start) > delta) {

    //}

}

void GeneratePassList(const CoordGeodetic& geo, const SGDP4& model, const Julian& date) {
    Observer obs(geo);
    Eci eci;
    model.FindPosition(eci, date);

    CoordTopographic topo = obs.GetLookAngle(eci);

    /*
     * set start and end date
     */
    Julian time0 = date;
    Julian time1 = date;
    time1.AddDay(10.0);

    /*
     * step throw period with 1 minute increments
     */
    for (Julian jd = date; jd <= time1; jd.AddMin(1.0)) {

    }
}

int main() {

    Tle tle = Tle("UK-DMC 2                ",
            "1 35683U 09041C   11089.11558659  .00000272  00000-0  54146-4 0  8712",
            "2 35683  98.0762 348.1067 0001434  99.8921 260.2456 14.69414094 89293");

    CoordGeodetic geo(Globals::Deg2Rad(51.360242), Globals::Deg2Rad(0.101473), 0.07);

    SGDP4 sgp4_model;
    sgp4_model.SetTle(tle);
    Julian date;

    GeneratePassList(geo, sgp4_model, date);


    return 0;
}

#if 0

http://olifantasia.com/projects/gnuradio/mdvh/weather_sat/weather_sat_scripts_without_capture_files_2010061701/decoding/poes-weather-hrpt-decoder/hrpt-decoder-1.0.0.2/satellite/predict/

xmnpda = 1440.0

    /* same formulas, but the one from predict is nicer */
    //sat->footprint = 2.0 * xkmper * acos (xkmper/sat->pos.w);
    sat->footprint = 12756.33 * acos (xkmper / (xkmper+sat->alt));
    age = sat->jul_utc - sat->jul_epoch;
    sat->orbit = (long) floor((sat->tle.xno * xmnpda/twopi +
                    age * sat->tle.bstar * ae) * age +
                    sat->tle.xmo/twopi) + sat->tle.revnum - 1;

bool TSat::IsGeostationary(void)
{
 /* This function returns a 1 if the satellite
    appears to be in a geostationary orbit

    Circular orbit at an altitude of 35 800 km over the equator.
    A satellite moving with the Earth's rotation in a geostationary
    orbit has a period of 23 hours, 56 minutes and 4 seconds.
 */
 double sma, aalt;

  if(meanmo == 0.0)
     return true;

  sma  = 331.25*exp(log(1440.0/meanmo)*(2.0/3.0));
  aalt = sma*(1.0+eccn)-xkmper;

  if(fabs(meanmo-omega_E) < 0.0005 || // allmost same speed as earth
     aalt > 35000)                    // altitude is over 35000 km
     return true;
  else
     return false;
}

// latitude in radians
bool TSat::DoesRise(double lat)
{
 /* This function returns a true if the satellite can ever rise
    above the horizon of the ground station.
*/
 double lin, sma, apogee;
 bool rc = false;

  if(meanmo == 0.0)
     return rc;
  else {
     lin = incl;

     if(lin >= 90.0)
        lin=180.0-lin;

     sma    = 331.25*exp(log(1440.0/meanmo)*(2.0/3.0));
     apogee = sma*(1.0+eccn)-xkmper;

     if((acos2(xkmper/(apogee+xkmper))+lin*deg2rad) > fabs(lat))
        rc = true;
     else
  	rc = false;
  }

 return rc;
}

#endif

