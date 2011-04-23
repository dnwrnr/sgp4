#include "CoordTopographic.h"

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




