#include "Observer.h"

#include "Globals.h"

/*
 * in degrees!
 */
Observer::Observer(const double latitude, const double longitude, const double altitude) {

    geo_.latitude = DegreesToRadians(latitude);
    geo_.longitude = DegreesToRadians(longitude);
    geo_.altitude = altitude;

    observers_eci_ = Eci(Julian(), geo_);
}

Observer::Observer(const CoordGeodetic &geo)
: geo_(geo) {

    observers_eci_ = Eci(Julian(), geo_);
}

Observer::~Observer(void) {
}

void Observer::UpdateObserversEci(const Julian &date) {

    /*
     * if date has changed, update for new date
     */
    if (observers_eci_.GetDate() != date) {
        observers_eci_ = Eci(date, geo_);
    }
}

/*
 * calculate lookangle between the observer and the passed in Eci object
 */
CoordTopographic Observer::GetLookAngle(const Eci &eci) {

    /*
     * update the observers Eci to match the time of the Eci passed in
     * if necessary
     */
    UpdateObserversEci(eci.GetDate());

    /*
     * calculate differences
     */
    Vector range_rate = eci.GetVelocity().Subtract(observers_eci_.GetVelocity());
    Vector range = eci.GetPosition().Subtract(observers_eci_.GetPosition());

    range.w = range.GetMagnitude();

    /*
     * Calculate Local Mean Sidereal Time for observers longitude
     */
    double theta = eci.GetDate().ToLocalMeanSiderealTime(geo_.longitude);

    double sin_lat = sin(geo_.latitude);
    double cos_lat = cos(geo_.latitude);
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    double top_s = sin_lat * cos_theta * range.x +
            sin_lat * sin_theta * range.y -
            cos_lat * range.z;
    double top_e = -sin_theta * range.x +
            cos_theta * range.y;
    double top_z = cos_lat * cos_theta * range.x +
            cos_lat * sin_theta * range.y +
            sin_lat * range.z;
    double az = atan(-top_e / top_s);

    if (top_s > 0.0)
        az += kPI;

    if (az < 0.0)
        az += 2.0 * kPI;

    double el = asin(top_z / range.w);
    double rate = range.Dot(range_rate) / range.w;

    /*
     * azimuth in radians
     * elevation in radians
     * range in km
     * range rate in km/s
     */
    CoordTopographic topo(az,
            el,
            range.w,
            rate);

    return topo;
}
