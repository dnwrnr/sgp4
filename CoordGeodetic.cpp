#include "CoordGeodetic.h"

CoordGeodetic::CoordGeodetic(const CoordGeodetic& b) {

    lat_ = b.lat_;
    lon_ = b.lon_;
    alt_ = b.alt_;
}

CoordGeodetic& CoordGeodetic::operator =(const CoordGeodetic& b) {

    if (this != &b) {
        lat_ = b.lat_;
        lon_ = b.lon_;
        alt_ = b.alt_;
    }

    return (*this);
}

bool CoordGeodetic::operator ==(const CoordGeodetic& b) const {

    if (lat_ == b.lat_ &&
            lon_ == b.lon_ &&
            alt_ == b.alt_) {
        return true;
    } else {
        return false;
    }
}

bool CoordGeodetic::operator !=(const CoordGeodetic& b) const {

    if (lat_ == b.lat_ &&
            lon_ == b.lon_ &&
            alt_ == b.alt_) {
        return false;
    } else {
        return true;
    }
}
