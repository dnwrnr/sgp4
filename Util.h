#ifndef UTIL_H_
#define UTIL_H_

#include "Globals.h"

#include <sstream>

namespace Util
{
    template
    <typename T>
    bool FromString(const std::string& str, T& val)
    {
        std::stringstream ss(str);
        return !(ss >> val).fail();
    }

    /*
     * always positive result
     * Mod(-3,4)= 1   fmod(-3,4)= -3
     */
    inline double Mod(const double x, const double y)
    {
        if (y == 0)
        {
            return x;
        }

        return x - y * floor(x / y);
    }

    inline double WrapNegPosPI(const double a)
    {
        return Mod(a + kPI, kTWOPI) - kPI;
    }
    
    inline double WrapTwoPI(const double a)
    {
        return Mod(a, kTWOPI);
    }

    inline double WrapNegPos180(const double a)
    {
        return Mod(a + 180.0, 360.0) - 180.0;
    }

    inline double Wrap360(const double a)
    {
        return Mod(a, 360.0);
    }

    inline double DegreesToRadians(const double degrees)
    {
        return degrees * kPI / 180.0;
    }

    inline double RadiansToDegrees(const double radians)
    {
        return radians * 180.0 / kPI;
    }

    inline double AcTan(const double sinx, const double cosx)
    {
        if (cosx == 0.0)
        {
            if (sinx > 0.0)
            {
                return kPI / 2.0;
            }
            else
            {
                return 3.0 * kPI / 2.0;
            }
        }
        else
        {
            if (cosx > 0.0)
            {
                return atan(sinx / cosx);
            }
            else
            {
                return kPI + atan(sinx / cosx);
            }
        }
    }

    void TrimLeft(std::string& s);
    void TrimRight(std::string& s);
    void Trim(std::string& s);
}

#endif
