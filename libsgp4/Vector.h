/*
 * Copyright 2013 Daniel Warner <contact@danrw.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#pragma once

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

namespace libsgp4
{

/**
 * @brief Generic vector
 *
 * Stores x, y, z, w
 */
struct Vector
{
public:

    /**
     * Default constructor
     */
    Vector() = default;

    /**
     * Constructor
     * @param argX x value
     * @param argY y value
     * @param argZ z value
     */
    Vector(double argX,
            double argY,
            double argZ)
        : x(argX), y(argY), z(argZ)
    {
    }

    /**
     * Constructor
     * @param argX x value
     * @param argY y value
     * @param argZ z value
     * @param argW w value
     */
    Vector(double argX,
            double argY,
            double argZ,
            double argW)
        : x(argX), y(argY), z(argZ), w(argW)
    {
    }

    /**
     * Subtract operator
     * @param v value to suctract from
     */
    Vector operator-(const Vector& v) const
    {
        return Vector(x - v.x,
                y - v.y,
                z - v.z,
                0.0);
    }

    /**
     * Calculates the magnitude of the vector
     * @returns magnitude of the vector
     */
    double Magnitude() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    /**
     * Calculates the dot product
     * @returns dot product
     */
    double Dot(const Vector& vec) const
    {
        return (x * vec.x) +
            (y * vec.y) +
            (z * vec.z);
    }

    /**
     * Converts this vector to a string
     * @returns this vector as a string
     */
    std::string ToString() const
    {
        std::stringstream ss;
        ss << std::right << std::fixed << std::setprecision(3);
        ss << "X: " << std::setw(9) << x;
        ss << ", Y: " << std::setw(9) << y;
        ss << ", Z: " << std::setw(9) << z;
        ss << ", W: " << std::setw(9) << w;
        return ss.str();
    }

    /** x value */
    double x{};
    /** y value */
    double y{};
    /** z value */
    double z{};
    /** w value */
    double w{};
};

inline std::ostream& operator<<(std::ostream& strm, const Vector& v)
{
    return strm << v.ToString();
}

} // namespace libsgp4
