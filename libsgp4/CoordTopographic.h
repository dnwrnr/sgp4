#ifndef COORDTOPOGRAPHIC_H_
#define COORDTOPOGRAPHIC_H_

#include "Util.h"

#include <string>
#include <sstream>
#include <iomanip>

/**
 * Stores a topographic position
 */
struct CoordTopographic
{
public:
    /**
     * Default constructor
     */
    CoordTopographic()
        : azimuth(0.0), 
        elevation(0.0),
        range(0.0),
        range_rate(0.0)
    {
    }

    /**
     * Constructor
     * @param[in] arg_azimuth azimuth in radians
     * @param[in] arg_elevation elevation in radians
     * @param[in] arg_range range in kilometers
     * @param[in] arg_range_rate range rate in kilometers per second
     */
    CoordTopographic(
            double arg_azimuth,
            double arg_elevation,
            double arg_range,
            double arg_range_rate)
        : azimuth(arg_azimuth),
        elevation(arg_elevation),
        range(arg_range),
        range_rate(arg_range_rate)
    {
    }

    /**
     * Copy constructor
     * @param[in] topo object to copy from
     */
    CoordTopographic(const CoordTopographic& topo)
    {
        azimuth = topo.azimuth;
        elevation = topo.elevation;
        range = topo.range;
        range_rate = topo.range_rate;
    }

    /**
     * Destructor
     */
    virtual ~CoordTopographic()
    {
    }

    /**
     * Assignment operator
     * @param[in] topo object to copy from
     */
    CoordTopographic& operator=(const CoordTopographic& topo)
    {
        if (this != &topo)
        {
            azimuth = topo.azimuth;
            elevation = topo.elevation;
            range = topo.range;
            range_rate = topo.range_rate;
        }
        return *this;
    }

    /**
     * Equality operator
     * @param[in] topo value to check
     * @returns whether the object is equal
     */
    bool operator==(const CoordTopographic& topo) const
    {
        return IsEqual(topo);
    }

    /**
     * Inequality operator
     * @param[in] topo the object to compare with
     * @returns whether the object is not equal
     */    
    bool operator !=(const CoordTopographic& topo) const
    {
        return !IsEqual(topo);
    }

    /**
     * Dump this object to a string
     * @returns string
     */
    std::string ToString() const
    {
        std::stringstream ss;
        ss << std::right << std::fixed << std::setprecision(3);
        ss << "Az: " << std::setw(8) << Util::RadiansToDegrees(azimuth);
        ss << ", El: " << std::setw(8) << Util::RadiansToDegrees(elevation);
        ss << ", Rng: " << std::setw(10) << range;
        ss << ", Rng Rt: " << std::setw(7) << range_rate;
        return ss.str();
    }

    /** azimuth in radians */
    double azimuth;
    /** elevations in radians */
    double elevation;
    /** range in kilometers */
    double range;
    /** range rate in kilometers per second */
    double range_rate;

private:
    bool IsEqual(const CoordTopographic& topo) const
    {
        bool equal = false;
        if (azimuth == topo.azimuth &&
                elevation == topo.elevation &&
                range == topo.range &&
                range_rate == topo.range_rate)
        {
            equal = true;
        }
        return equal;
    }
};


inline std::ostream& operator<<(std::ostream& strm, const CoordTopographic& t)
{
    return strm << t.ToString();
}

#endif
