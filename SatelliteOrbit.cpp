#include "SatelliteOrbit.h"

SatelliteOrbit::SatelliteOrbit(void) {
}

SatelliteOrbit::~SatelliteOrbit(void) {
}

void SatelliteOrbit::SetTle(const Tle& tle) {

    sgdp4_.SetTle(tle);
}

bool SatelliteOrbit::IsGeostationary() {
#if 0
    if (sgdp4_.MeanMotion() == 0.0)
        return true;

    /*
    radius of apogee 
    the distance from the centre of the planet to the point in the orbit furthest away from the planet
    */
    const double apogee_altitude = sgdp4_.RecoveredSemiMajorAxis() * (1.0 + sgdp4_.Eccentricity()) - Globals::XKMPER();

    /*
     * check if almost same speed as earth
     * or altitude is over 35000 km
     */
    if (fabs(sgdp4_.MeanMotion() - Globals::OMEGA_E()) < 0.0005 || apogee_altitude > 35000)
        return true;
    else
#endif
        return false;
}

unsigned int SatelliteOrbit::GetOrbitNumber(const Julian& jul) const{
#if 0
    double diff = jul.SpanMin(sgdp4_.Epoch());
    
    return (unsigned int)floor((sgdp4_.MeanMotion() * 1440.0 / Globals::TWOPI() +
        diff * sgdp4_.BStar() * Globals::AE()) * diff +
        sgdp4_.MeanAnomoly() / Globals::TWOPI()) + sgdp4_.OrbitNumber() - 1.0;
#endif
    return 0;
}

#if 0
/* same formulas, but the one from predict is nicer */
//sat->footprint = 2.0 * xkmper * acos (xkmper/sat->pos.w);
sat->footprint = 12756.33 * acos(xkmper / (xkmper + sat->alt));
age = sat->jul_utc - sat->jul_epoch;
sat->orbit = (long) floor((sat->tle.xno * xmnpda / twopi +
        age * sat->tle.bstar * ae) * age +
        sat->tle.xmo / twopi) + sat->tle.revnum - 1;
#endif