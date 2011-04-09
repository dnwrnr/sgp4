#include "Coord.h"

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

CoordTopographic::CoordTopographic(const CoordTopographic& b) {

    azimuth_ = b.azimuth_;
    elevation_ = b.elevation_;
    range_ = b.range_;
    range_rate_ = b.range_rate_;
}

CoordTopographic& CoordTopographic::operator =(const CoordTopographic& b) {

    if (this != &b) {
        azimuth_ = b.azimuth_;
        elevation_ = b.elevation_;
        range_ = b.range_;
        range_rate_ = b.range_rate_;
    }

    return (*this);
}

bool CoordTopographic::operator ==(const CoordTopographic& b) const {

    if (azimuth_ == b.azimuth_ &&
            elevation_ == b.elevation_ &&
            range_ == b.range_ &&
            range_rate_ == b.range_rate_) {
        return true;
    } else {
        return false;
    }
}

bool CoordTopographic::operator !=(const CoordTopographic& b) const {

    if (azimuth_ == b.azimuth_ &&
            elevation_ == b.elevation_ &&
            range_ == b.range_ &&
            range_rate_ == b.range_rate_) {
        return false;
    } else {
        return true;
    }
}




