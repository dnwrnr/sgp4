#ifndef COORDGEODETIC_H_
#define COORDGEODETIC_H_

#include "Util.h"

#include <string>
#include <sstream>
#include <iomanip>

/**
 * Stores a geodetic position
 */
struct CoordGeodetic
{
public:
    /**
     * Default constructor
     */
    CoordGeodetic()
        : latitude(0.0),
        longitude(0.0),
        altitude(0.0)
    {
    }

    /**
     * Constructor
     * @param[in] arg_latitude the latitude in degrees
     * @param[in] arg_longitude the longitude in degrees
     * @param[in] arg_altitude the altitude in kilometers
     */
    CoordGeodetic(
            double arg_latitude,
            double arg_longitude,
            double arg_altitude)
    {
        latitude = Util::DegreesToRadians(arg_latitude);
        longitude = Util::DegreesToRadians(arg_longitude);
        altitude = arg_altitude;
    }

    /**
     * Copy constructor
     * @param[in] geo object to copy from
     */
    CoordGeodetic(const CoordGeodetic& geo)
    {
        latitude = geo.latitude;
        longitude = geo.longitude;
        altitude = geo.altitude;
    }

    /**
     * Destructor
     */
    virtual ~CoordGeodetic()
    {
    }

    /**
     * Assignment operator
     * @param[in] geo object to copy from
     */
    CoordGeodetic& operator=(const CoordGeodetic& geo)
    {
        if (this != &geo)
        {
            latitude = geo.latitude;
            longitude = geo.longitude;
            altitude = geo.altitude;
        }
        return *this;
    }

    /**
     * Equality operator
     * @param[in] geo the object to compare with
     * @returns whether the object is equal
     */
    bool operator==(const CoordGeodetic& geo) const
    {
        return IsEqual(geo);
    }

    /**
     * Inequality operator
     * @param[in] geo the object to compare with
     * @returns whether the object is not equal
     */
    bool operator!=(const CoordGeodetic& geo) const
    {
        return !IsEqual(geo);
    }

    /**
     * Dump this object to a string
     * @returns string
     */
    std::string ToString() const
    {
        std::stringstream ss;
        ss << std::right << std::fixed << std::setprecision(3);
        ss << "Lat: " << std::setw(7) << Util::RadiansToDegrees(latitude);
        ss << ", Lon: " << std::setw(7) << Util::RadiansToDegrees(longitude);
        ss << ", Alt: " << std::setw(9) << altitude;
        return ss.str();
    }

    /** latitude in radians (-PI >= latitude < PI) */
    double latitude;
    /** latitude in radians (-PI/2 >= latitude <= PI/2) */
    double longitude;
    /** altitude in kilometers */
    double altitude;

private:
    bool IsEqual(const CoordGeodetic& geo) const
    {
        bool equal = false;
        if (latitude == geo.latitude &&
                longitude == geo.longitude &&
                altitude == geo.altitude)
        {
            equal = false;
        }
        return equal;
    }
};

inline std::ostream& operator<<(std::ostream& strm, const CoordGeodetic& g)
{
    return strm << g.ToString();
}

#endif
