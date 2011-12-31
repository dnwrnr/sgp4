#ifndef COORDGEODETIC_H_
#define COORDGEODETIC_H_

#include "Globals.h"
#include "Util.h"

#include <string>
#include <sstream>
#include <iomanip>

struct CoordGeodetic
{
public:
    CoordGeodetic()
        : latitude(0.0), longitude(0.0), altitude(0.0)
    {
    }

    /*
     * default is in degrees 
     */
    CoordGeodetic(double lat, double lon, double alt, bool radians = false)
    {
        if (radians)
        {
            latitude = lat;
            longitude = lon;
        }
        else
        {
            latitude = Util::DegreesToRadians(lat);
            longitude = Util::DegreesToRadians(lon);
        }
        altitude = alt;
    }

    CoordGeodetic(const CoordGeodetic& g)
    {
        latitude = g.latitude;
        longitude = g.longitude;
        altitude = g.altitude;
    }

    virtual ~CoordGeodetic()
    {
    }

    CoordGeodetic& operator=(const CoordGeodetic& g)
    {
        if (this != &g)
        {
            latitude = g.latitude;
            longitude = g.longitude;
            altitude = g.altitude;
        }
        return *this;
    }

    bool operator==(const CoordGeodetic& g) const
    {
        return IsEqual(g);
    }

    bool operator!=(const CoordGeodetic& g) const
    {
        return !IsEqual(g);
    }

    std::string ToString() const
    {
        std::stringstream ss;
        ss << std::right << std::fixed << std::setprecision(2);
        ss << "Lat: " << std::setw(7) << Util::RadiansToDegrees(latitude);
        ss << ", Lon: " << std::setw(7) << Util::RadiansToDegrees(longitude);
        ss << ", Alt: " << std::setw(9) << altitude;
        return ss.str();
    }

    /*
     * radians (north positive, south negative)
     */
    double latitude;
    /*
     * radians (east positive, west negative)
     */
    double longitude;
    /*
     * kilometers
     */
    double altitude;

protected:
    bool IsEqual(const CoordGeodetic& g) const
    {
        if (latitude == g.latitude && longitude == g.longitude &&
            altitude == g.altitude)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
};

inline std::ostream& operator<<(std::ostream& strm, const CoordGeodetic& g)
{
    return strm << g.ToString();
}

#endif

