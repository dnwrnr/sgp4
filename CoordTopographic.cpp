#include "CoordTopographic.h"

CoordTopographic::CoordTopographic(const CoordTopographic& b) {

    azimuth = b.azimuth;
    elevation = b.elevation;
    range = b.range;
    range_rate = b.range_rate;
}

CoordTopographic& CoordTopographic::operator =(const CoordTopographic& b) {

    if (this != &b) {
        azimuth = b.azimuth;
        elevation = b.elevation;
        range = b.range;
        range_rate = b.range_rate;
    }

    return (*this);
}

bool CoordTopographic::operator ==(const CoordTopographic& b) const {

    if (azimuth == b.azimuth &&
            elevation == b.elevation &&
            range == b.range &&
            range_rate == b.range_rate) {
        return true;
    } else {
        return false;
    }
}

bool CoordTopographic::operator !=(const CoordTopographic& b) const {

    if (azimuth == b.azimuth &&
            elevation == b.elevation &&
            range == b.range &&
            range_rate == b.range_rate) {
        return false;
    } else {
        return true;
    }
}




