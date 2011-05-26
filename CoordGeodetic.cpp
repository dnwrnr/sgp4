#include "CoordGeodetic.h"

CoordGeodetic::CoordGeodetic(const CoordGeodetic& b) {

    latitude = b.latitude;
    longitude = b.longitude;
    altitude = b.altitude;
}

CoordGeodetic& CoordGeodetic::operator =(const CoordGeodetic& b) {

    if (this != &b) {
        latitude = b.latitude;
        longitude = b.longitude;
        altitude = b.altitude;
    }

    return (*this);
}

bool CoordGeodetic::operator ==(const CoordGeodetic& b) const {

    if (latitude == b.latitude &&
            longitude == b.longitude &&
            altitude == b.altitude) {
        return true;
    } else {
        return false;
    }
}

bool CoordGeodetic::operator !=(const CoordGeodetic& b) const {

    if (latitude == b.latitude &&
            longitude == b.longitude &&
            altitude == b.altitude) {
        return false;
    } else {
        return true;
    }
}
