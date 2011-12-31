#include "Eci.h"

#include "Util.h"

void Eci::ToEci(const Julian& date, const CoordGeodetic &g)
{
    /*
     * set date
     */
    date_ = date;

    static const double mfactor = kTWOPI * (kOMEGA_E / kSECONDS_PER_DAY);

    /*
     * Calculate Local Mean Sidereal Time for observers longitude
     */
    const double theta = date_.ToLocalMeanSiderealTime(g.longitude);

    /*
     * take into account earth flattening
     */
    const double c = 1.0
        / sqrt(1.0 + kF * (kF - 2.0) * pow(sin(g.latitude), 2.0));
    const double s = pow(1.0 - kF, 2.0) * c;
    const double achcp = (kXKMPER * c + g.altitude) * cos(g.latitude);

    /*
     * X position in km
     * Y position in km
     * Z position in km
     * W magnitude in km
     */
    position_.x = achcp * cos(theta);
    position_.y = achcp * sin(theta);
    position_.z = (kXKMPER * s + g.altitude) * sin(g.latitude);
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

CoordGeodetic Eci::ToGeodetic() const
{
    const double theta = Util::AcTan(position_.y, position_.x);

    // 0 >= lon < 360
    // const double lon = Fmod2p(theta - date_.ToGreenwichSiderealTime());
    // 180 >= lon < 180
    const double lon = Util::WrapNegPosPI(theta - date_.ToGreenwichSiderealTime());

    const double r = sqrt((position_.x * position_.x)
            + (position_.y * position_.y));
    
    static const double e2 = kF * (2.0 - kF);

    double lat = Util::AcTan(position_.z, r);
    double phi = 0.0;
    double c = 0.0;
    int cnt = 0;

    do
    {
        phi = lat;
        const double sinphi = sin(phi);
        c = 1.0 / sqrt(1.0 - e2 * sinphi * sinphi);
        lat = Util::AcTan(position_.z + kXKMPER * c * e2 * sinphi, r);
        cnt++;
    }
    while (fabs(lat - phi) >= 1e-10 && cnt < 10);

    const double alt = r / cos(lat) - kXKMPER * c;

    return CoordGeodetic(lat, lon, alt, true);
}
