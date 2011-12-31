#ifndef COORDTOPOGRAPHIC_H_
#define COORDTOPOGRAPHIC_H_

#include "Globals.h"
#include "Util.h"

#include <iostream>
#include <sstream>
#include <iomanip>

struct CoordTopographic
{
public:
    CoordTopographic()
    : azimuth(0.0), elevation(0.0), range(0.0), range_rate(0.0)
    {
    }

    CoordTopographic(double az, double el, double rnge, double rnge_rate)
    : azimuth(az), elevation(el), range(rnge), range_rate(rnge_rate)
    {
    }

    CoordTopographic(const CoordTopographic& t)
    {
        azimuth = t.azimuth;
        elevation = t.elevation;
        range = t.range;
        range_rate = t.range_rate;
    }

    virtual ~CoordTopographic()
    {
    };

    CoordTopographic& operator=(const CoordTopographic& t)
    {
        if (this != &t) {
            azimuth = t.azimuth;
            elevation = t.elevation;
            range = t.range;
            range_rate = t.range_rate;
        }
        return *this;
    }

    bool operator==(const CoordTopographic& t) const
    {
        return IsEqual(t);
    }

    bool operator !=(const CoordTopographic& t) const
    {
        return !IsEqual(t);
    }

    std::string ToString() const
    {
        std::stringstream ss;
        ss << std::right << std::fixed << std::setprecision(2);
        ss << "Az: " << std::setw(7) << Util::RadiansToDegrees(azimuth);
        ss << ", El: " << std::setw(7) << Util::RadiansToDegrees(elevation);
        ss << ", Rng: " << std::setw(9) << range;
        ss << ", Rng Rt: " << std::setw(6) << range_rate;
        return ss.str();
    }

    /*
     * radians
     */
    double azimuth;
    /*
     * radians
     */
    double elevation;
    /*
     * kilometers
     */
    double range;
    /*
     * kilometers / second
     */
    double range_rate;

protected:
    bool IsEqual(const CoordTopographic& t) const
    {
       if (azimuth == t.azimuth && elevation == t.elevation &&
                range == t.range && range_rate == t.range_rate)
       {
            return true;
       }
       else
       {
           return false;
       }
    }
};


inline std::ostream& operator<<(std::ostream& strm, const CoordTopographic& t)
{
    return strm << t.ToString();
}

#endif

