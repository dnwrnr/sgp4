#include "Vector.h"

double Vector::GetMagnitude() const {
    return sqrt(x * x + y * y + z * z);
}

/*
 * subtract (this) - (v)
 * and return result
 */
Vector Vector::Subtract(const Vector& vec) const {
    return Vector(x - vec.x,
            y - vec.y,
            z - vec.z,
            0.0);
}

double Vector::Dot(const Vector& vec) const {

    return (x * vec.x) +
            (y * vec.y) +
            (z * vec.z);
}
