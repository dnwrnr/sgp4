#include "Eci.h"

Eci::Eci(const Julian &date, const CoordGeodetic &geo)
: date_(date) {

    static const double mfactor = Globals::TWOPI() * (Globals::OMEGA_E() / Globals::SEC_PER_DAY());
    const double latitude = geo.GetLatitude();
    const double longitude = geo.GetLongitude();
    const double altitude = geo.GetAltitude();

    /*
     * Calculate Local Mean Sidereal Time for observers longitude
     */
    const double theta = date_.ToLocalMeanSiderealTime(longitude);

    /*
     * take into account earth flattening
     */
    const double c = 1.0 / sqrt(1.0 + Globals::F() * (Globals::F() - 2.0) * pow(sin(latitude), 2.0));
    const double s = pow(1.0 - Globals::F(), 2.0) * c;
    const double achcp = (Globals::XKMPER() * c + altitude) * cos(latitude);

    /*
     * X position in km
     * Y position in km
     * Z position in km
     * W magnitude in km
     */
    position_.SetX(achcp * cos(theta));
    position_.SetY(achcp * sin(theta));
    position_.SetZ((Globals::XKMPER() * s + altitude) * sin(latitude));
    position_.SetW(position_.GetMagnitude());

    /*
     * X velocity in km/s
     * Y velocity in km/s
     * Z velocity in km/s
     * W magnitude in km/s
     */
    velocity_.SetX(-mfactor * position_.GetY());
    velocity_.SetY(mfactor * position_.GetX());
    velocity_.SetZ(0.0);
    velocity_.SetW(velocity_.GetMagnitude());
}

Eci::Eci(const Julian &date, const Vector &position, const Vector &velocity)
: date_(date), position_(position), velocity_(velocity) {

}

Eci::~Eci(void) {
}

CoordGeodetic Eci::ToGeodetic() const {

    const double theta = Globals::AcTan(position_.GetY(), position_.GetX());
    /*
     * changes lon to 0>= and <360
     * const double lon = Globals::Fmod2p(theta - date_.ToGreenwichSiderealTime());
     */
    const double lon = fmod(theta - date_.ToGreenwichSiderealTime(), Globals::TWOPI());
    const double r = sqrt((position_.GetX() * position_.GetX()) + (position_.GetY() * position_.GetY()));
    static const double e2 = Globals::F() * (2.0 - Globals::F());

    double lat = Globals::AcTan(position_.GetZ(), r);
    double phi = 0.0;
    double c = 0.0;
    int cnt = 0;

    do {
        phi = lat;
        const double sinphi = sin(phi);
        c = 1.0 / sqrt(1.0 - e2 * sinphi * sinphi);
        lat = Globals::AcTan(position_.GetZ() + Globals::XKMPER() * c * e2 * sinphi, r);
        cnt++;
    } while (fabs(lat - phi) >= 1e-10 && cnt < 10);

    const double alt = r / cos(lat) - Globals::XKMPER() * c;

    return CoordGeodetic(lat, lon, alt);
}
