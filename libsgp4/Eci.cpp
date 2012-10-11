#include "Eci.h"

#include "Globals.h"
#include "Util.h"

/**
 * Converts a DateTime and Geodetic position to Eci coordinates
 * @param[in] dt the date
 * @param[in] geo the geodetic position
 */
void Eci::ToEci(const DateTime& dt, const CoordGeodetic &geo)
{
    /*
     * set date
     */
    m_dt = dt;

    static const double mfactor = kTWOPI * (kOMEGA_E / kSECONDS_PER_DAY);

    /*
     * Calculate Local Mean Sidereal Time for observers longitude
     */
    const double theta = m_dt.ToLocalMeanSiderealTime(geo.longitude);

    /*
     * take into account earth flattening
     */
    const double c = 1.0
        / sqrt(1.0 + kF * (kF - 2.0) * pow(sin(geo.latitude), 2.0));
    const double s = pow(1.0 - kF, 2.0) * c;
    const double achcp = (kXKMPER * c + geo.altitude) * cos(geo.latitude);

    /*
     * X position in km
     * Y position in km
     * Z position in km
     * W magnitude in km
     */
    m_position.x = achcp * cos(theta);
    m_position.y = achcp * sin(theta);
    m_position.z = (kXKMPER * s + geo.altitude) * sin(geo.latitude);
    m_position.w = m_position.GetMagnitude();

    /*
     * X velocity in km/s
     * Y velocity in km/s
     * Z velocity in km/s
     * W magnitude in km/s
     */
    m_velocity.x = -mfactor * m_position.y;
    m_velocity.y = mfactor * m_position.x;
    m_velocity.z = 0.0;
    m_velocity.w = m_velocity.GetMagnitude();
}

/**
 * @returns the position in geodetic form
 */
CoordGeodetic Eci::ToGeodetic() const
{
    const double theta = Util::AcTan(m_position.y, m_position.x);

    const double lon = Util::WrapNegPosPI(theta
            - m_dt.ToGreenwichSiderealTime());

    const double r = sqrt((m_position.x * m_position.x)
            + (m_position.y * m_position.y));
    
    static const double e2 = kF * (2.0 - kF);

    double lat = Util::AcTan(m_position.z, r);
    double phi = 0.0;
    double c = 0.0;
    int cnt = 0;

    do
    {
        phi = lat;
        const double sinphi = sin(phi);
        c = 1.0 / sqrt(1.0 - e2 * sinphi * sinphi);
        lat = Util::AcTan(m_position.z + kXKMPER * c * e2 * sinphi, r);
        cnt++;
    }
    while (fabs(lat - phi) >= 1e-10 && cnt < 10);

    const double alt = r / cos(lat) - kXKMPER * c;

    return CoordGeodetic(lat, lon, alt);
}
