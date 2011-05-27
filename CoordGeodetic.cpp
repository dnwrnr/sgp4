#include "CoordGeodetic.h"

#include <sstream>
#include <iomanip>

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

std::ostream& operator<< (std::ostream& stream, const CoordGeodetic& geo) {
    std::stringstream out;
    out << std::right << std::fixed << std::setprecision(2);
    out << "Lat: " << std::setw(7) << geo.latitude;
    out << ", Lon: " << std::setw(7) << geo.longitude;
    out << ", Alt: " << std::setw(9) << geo.altitude;
    stream << out.str();
    return stream;
}
