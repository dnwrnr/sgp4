#ifndef VECTOR_H_
#define VECTOR_H_

#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

struct Vector
{
public:

    Vector()
        : x(0.0), y(0.0), z(0.0), w(0.0)
    {
    }

    Vector(double xx, double yy, double zz)
        : x(xx), y(yy), z(zz), w(0.0)
    {
    }

    Vector(double xx, double yy, double zz, double ww)
        : x(xx), y(yy), z(zz), w(ww)
    {
    }
    
    Vector(const Vector& v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
        w = v.w;
    }

    virtual ~Vector()
    {
    }

    Vector& operator=(const Vector& v)
    {
        if (this != &v)
        {
            x = v.x;
            y = v.y;
            z = v.z;
            w = v.w;
        }
        return *this;
    }

    double GetMagnitude() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    Vector Subtract(const Vector& v) const
    {
        /*
         * subtract (this) - (v)
         * and return result
         */
        return Vector(x - v.x,
                y - v.y,
                z - v.z,
                0.0);
    }

    double Dot(const Vector& vec) const
    {
        return (x * vec.x) +
            (y * vec.y) +
            (z * vec.z);
    }

    std::string ToString() const
    {
        std::stringstream ss;
        ss << std::right << std::fixed << std::setprecision(2);
        ss << "X: " << std::setw(8) << x;
        ss << ", Y: " << std::setw(8) << y;
        ss << ", Z: " << std::setw(8) << z;
        ss << ", W: " << std::setw(8) << w;
        return ss.str();
    }

    double x;
    double y;
    double z;
    double w;
};

inline std::ostream& operator<<(std::ostream& strm, const Vector& v)
{
    return strm << v.ToString();
}

#endif
