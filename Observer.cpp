#include "Observer.h"

#include "Globals.h"

/*
 * in degrees!
 */
Observer::Observer(const double latitude, const double longitude, const double altitude) {
    geo_.SetLatitude(Globals::Deg2Rad(latitude));
    geo_.SetLongitude(Globals::Deg2Rad(longitude));
    geo_.SetAltitude(altitude);

    observers_eci_ = Eci(Julian(), geo_);
}

Observer::Observer(const CoordGeodetic &geo)
: geo_(geo), observers_eci_(Julian(), geo) {

}

Observer::~Observer(void) {
}

void Observer::UpdateObserversEci(const Julian &date) {

    if (observers_eci_.GetDate() != date) {
        observers_eci_ = Eci(date, geo_);
    }
}

CoordTopographic Observer::GetLookAngle(const Eci &eci) {

    Julian date = eci.GetDate();

    /*
     * update observers eci value if the date is different to the eci passed in
     */
    UpdateObserversEci(date);

    Vector range_rate(eci.GetVelocity().GetX() - observers_eci_.GetVelocity().GetX(),
            eci.GetVelocity().GetY() - observers_eci_.GetVelocity().GetY(),
            eci.GetVelocity().GetZ() - observers_eci_.GetVelocity().GetZ());

    const double x = eci.GetPosition().GetX() - observers_eci_.GetPosition().GetX();
    const double y = eci.GetPosition().GetY() - observers_eci_.GetPosition().GetY();
    const double z = eci.GetPosition().GetZ() - observers_eci_.GetPosition().GetZ();
    const double w = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));

    Vector range(x, y, z, w);

    // The site's Local Mean Sidereal Time at the time of interest.
    double theta = date.ToLocalMeanSiderealTime(geo_.GetLongitude());

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
    double rate = (range.GetX() * range_rate.GetX() +
            range.GetY() * range_rate.GetY() +
            range.GetZ() * range_rate.GetZ()) / range.GetW();

    CoordTopographic topo(az, // azimuth,   radians
            el, // elevation, radians
            range.GetW(), // range, km
            rate); // rate,  km / sec

    return topo;
}
