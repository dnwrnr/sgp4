#include "Observer.h"

#include "Globals.h"

/*
 * in degrees!
 */
Observer::Observer(const double latitude, const double longitude, const double altitude) {

    geo_.SetLatitude(Globals::Deg2Rad(latitude));
    geo_.SetLongitude(Globals::Deg2Rad(longitude));
    geo_.SetAltitude(altitude);

    UpdateObserversEci(Julian());
}

Observer::Observer(const CoordGeodetic &geo)
: geo_(geo) {

    UpdateObserversEci(Julian());
}

Observer::~Observer(void) {
}

void Observer::UpdateObserversEci(const Julian &date) {

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
     */
    UpdateObserversEci(eci.GetDate());

    /*
     * calculate differences
     */
    Vector range_rate = eci.GetVelocity().Subtract(observers_eci_.GetVelocity());
    Vector range = eci.GetPosition().Subtract(observers_eci_.GetPosition());

    range.SetW(range.GetMagnitude());

    /*
     * Calculate Local Mean Sidereal Time for observers longitude
     */
    double theta = eci.GetDate().ToLocalMeanSiderealTime(geo_.GetLongitude());

    double sin_lat = sin(geo_.GetLatitude());
    double cos_lat = cos(geo_.GetLatitude());
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    double top_s = sin_lat * cos_theta * range.GetX() +
            sin_lat * sin_theta * range.GetY() -
            cos_lat * range.GetZ();
    double top_e = -sin_theta * range.GetX() +
            cos_theta * range.GetY();
    double top_z = cos_lat * cos_theta * range.GetX() +
            cos_lat * sin_theta * range.GetY() +
            sin_lat * range.GetZ();
    double az = atan(-top_e / top_s);

    if (top_s > 0.0)
        az += Globals::PI();

    if (az < 0.0)
        az += 2.0 * Globals::PI();

    double el = asin(top_z / range.GetW());
    double rate = range.Dot(range_rate) / range.GetW();

    /*
     * azimuth in radians
     * elevation in radians
     * range in km
     * range rate in km/s
     */
    CoordTopographic topo(az,
            el,
            range.GetW(),
            rate);

    return topo;
}
