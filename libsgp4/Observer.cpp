#include "Observer.h"

#include "CoordTopocentric.h"

/*
 * calculate lookangle between the observer and the passed in Eci object
 */
CoordTopocentric Observer::GetLookAngle(const Eci &eci)
{
    /*
     * update the observers Eci to match the time of the Eci passed in
     * if necessary
     */
    Update(eci.GetDateTime());

    /*
     * calculate differences
     */
    Vector range_rate = eci.Velocity() - m_eci.Velocity();
    Vector range = eci.Position() - m_eci.Position();

    range.w = range.Magnitude();

    /*
     * Calculate Local Mean Sidereal Time for observers longitude
     */
    double theta = eci.GetDateTime().ToLocalMeanSiderealTime(m_geo.longitude);

    double sin_lat = sin(m_geo.latitude);
    double cos_lat = cos(m_geo.latitude);
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    double top_s = sin_lat * cos_theta * range.x
        + sin_lat * sin_theta * range.y - cos_lat * range.z;
    double top_e = -sin_theta * range.x
        + cos_theta * range.y;
    double top_z = cos_lat * cos_theta * range.x 
        + cos_lat * sin_theta * range.y + sin_lat * range.z;
    double az = atan(-top_e / top_s);

    if (top_s > 0.0)
    {
        az += kPI;
    }

    if (az < 0.0)
    {
        az += 2.0 * kPI;
    }

    double el = asin(top_z / range.w);
    double rate = range.Dot(range_rate) / range.w;

    /*
     * azimuth in radians
     * elevation in radians
     * range in km
     * range rate in km/s
     */
    return CoordTopocentric(az,
            el,
            range.w,
            rate);
}
