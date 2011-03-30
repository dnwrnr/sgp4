#include "Observer.h"

#include "Globals.h"

/*
 * in degrees!
 */
Observer::Observer(const double latitude, const double longitude, const double altitude) {
    geo_.SetLatitude(Globals::Deg2Rad(latitude));
    geo_.SetLongitude(Globals::Deg2Rad(longitude));
    geo_.SetAltitude(altitude);
}

Observer::Observer(const CoordGeodetic &geo)
: geo_(geo) {

}

Observer::~Observer(void) {
}
