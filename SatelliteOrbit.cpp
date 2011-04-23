#include "SatelliteOrbit.h"

SatelliteOrbit::SatelliteOrbit(void) {
}

SatelliteOrbit::~SatelliteOrbit(void) {
}

void SatelliteOrbit::SetTle(const Tle& tle) {

    sgp4_.SetTle(tle);
}

bool SatelliteOrbit::IsGeostationary() {

    if (sgp4_.MeanMotion() == 0.0)
        return true;

    /*
    radius of apogee
    the distance from the centre of the planet to the point in the orbit furthest away from the planet
     */
    const double apogee_altitude = sgp4_.RecoveredSemiMajorAxis() * (1.0 + sgp4_.Eccentricity()) - Globals::XKMPER();

    /*
     * check if almost same speed as earth
     * or altitude is over 35000 km
     */
    if (fabs(sgp4_.MeanMotion() - Globals::OMEGA_E()) < 0.0005 || apogee_altitude > 35000)
        return true;
    else
        return false;
}

unsigned int SatelliteOrbit::GetOrbitNumber(const Julian& jul) const {

    double diff = jul.SpanMin(sgp4_.Epoch());

    return static_cast<unsigned int> (floor((sgp4_.MeanMotion() * 1440.0 / Globals::TWOPI() +
            diff * sgp4_.BStar() * Globals::AE()) * diff +
            sgp4_.MeanAnomoly() / Globals::TWOPI())) + sgp4_.OrbitNumber() - 1;
}

double SatelliteOrbit::Footprint(const double& altitude) {

    if (altitude > 0)
        return 2.0 * Globals::XKMPER() * acos(Globals::XKMPER() / (Globals::XKMPER() + altitude));
    else
        return 0.0;
}