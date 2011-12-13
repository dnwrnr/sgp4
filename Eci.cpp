#include "Eci.h"

#include "Globals.h"

Eci::Eci(const Julian &date, const CoordGeodetic &geo)
: date_(date) {

    static const double mfactor = kTWOPI * (kOMEGA_E / kSECONDS_PER_DAY);
    const double latitude = geo.latitude;
    const double longitude = geo.longitude;
    const double altitude = geo.altitude;

    /*
     * Calculate Local Mean Sidereal Time for observers longitude
     */
    const double theta = date_.ToLocalMeanSiderealTime(longitude);

    /*
     * take into account earth flattening
     */
    const double c = 1.0 / sqrt(1.0 + kF * (kF - 2.0) * pow(sin(latitude), 2.0));
    const double s = pow(1.0 - kF, 2.0) * c;
    const double achcp = (kXKMPER * c + altitude) * cos(latitude);

    /*
     * X position in km
     * Y position in km
     * Z position in km
     * W magnitude in km
     */
    position_.x = achcp * cos(theta);
    position_.y = achcp * sin(theta);
    position_.z = (kXKMPER * s + altitude) * sin(latitude);
    position_.w = position_.GetMagnitude();

    /*
     * X velocity in km/s
     * Y velocity in km/s
     * Z velocity in km/s
     * W magnitude in km/s
     */
    velocity_.x = -mfactor * position_.y;
    velocity_.y = mfactor * position_.x;
    velocity_.z = 0.0;
    velocity_.w = velocity_.GetMagnitude();
}

Eci::Eci(const Julian &date, const Vector &position)
: date_(date), position_(position) {

}

Eci::Eci(const Julian &date, const Vector &position, const Vector &velocity)
: date_(date), position_(position), velocity_(velocity) {

}

Eci::~Eci(void) {
}

CoordGeodetic Eci::ToGeodetic() const {

    const double theta = AcTan(position_.y, position_.x);

    // 0 >= lon < 360
    // const double lon = Fmod2p(theta - date_.ToGreenwichSiderealTime());
    // 180 >= lon < 180
    const double lon = fmod(theta - date_.ToGreenwichSiderealTime(), kPI);

    const double r = sqrt((position_.x * position_.x) + (position_.y * position_.y));
    static const double e2 = kF * (2.0 - kF);

    double lat = AcTan(position_.z, r);
    double phi = 0.0;
    double c = 0.0;
    int cnt = 0;

    do {
        phi = lat;
        const double sinphi = sin(phi);
        c = 1.0 / sqrt(1.0 - e2 * sinphi * sinphi);
        lat = AcTan(position_.z + kXKMPER * c * e2 * sinphi, r);
        cnt++;
    } while (fabs(lat - phi) >= 1e-10 && cnt < 10);

    const double alt = r / cos(lat) - kXKMPER * c;

    return CoordGeodetic(lat, lon, alt);
}
