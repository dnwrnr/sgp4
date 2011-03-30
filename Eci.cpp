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
    const double theta = date_.ToLMST(longitude);

    const double c = 1.0 / sqrt(1.0 + Globals::F() * (Globals::F() - 2.0) * pow(sin(latitude), 2.0));
    const double s = pow(1.0 - Globals::F(), 2.0) * c;
    const double achcp = (Globals::XKMPER() * c + altitude) * cos(latitude);

    position_.SetX(achcp * cos(theta));
    position_.SetY(achcp * sin(theta));

    position_.SetZ((Globals::XKMPER() * s + altitude) * sin(latitude));
    position_.SetW(sqrt(pow(position_.GetX(), 2.0) + pow(position_.GetY(), 2.0) + pow(position_.GetZ(), 2.0)));

    velocity_.SetX(-mfactor * position_.GetY());
    velocity_.SetY(mfactor * position_.GetX());
    velocity_.SetZ(0.0);
    velocity_.SetW(sqrt(pow(velocity_.GetX(), 2.0) + pow(velocity_.GetY(), 2.0)));
}

Eci::Eci(const Julian &date, const Vector &position, const Vector &velocity)
: date_(date), position_(position), velocity_(velocity) {

}

Eci::~Eci(void) {
}
