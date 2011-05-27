#include "CoordTopographic.h"

#include "Globals.h"

#include <sstream>
#include <iomanip>

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

std::ostream& operator<< (std::ostream& stream, const CoordTopographic& topo) {
    std::stringstream out;
    out << std::right << std::fixed << std::setprecision(2);
    out << "Az: " << std::setw(7) << RadiansToDegrees(topo.azimuth);
    out << ", El: " << std::setw(7) << RadiansToDegrees(topo.elevation);
    out << ", Range: " << std::setw(9) << topo.range;
    out << ", RangeRate: " << std::setw(6) << topo.range_rate;
    stream << out.str();
    return stream;
}


