#include "Vector.h"

double Vector::GetMagnitude() const {
    return sqrt(x_ * x_ + y_ * y_ + z_ * z_);
}

/*
 * subtract (this) - (v)
 * and return result
 */
Vector Vector::Subtract(const Vector& vec) const {
    return Vector(x_ - vec.x_,
            y_ - vec.y_,
            z_ - vec.z_,
            0.0);
}

double Vector::Dot(const Vector& vec) const {

    return (x_ * vec.x_) +
            (y_ * vec.y_) +
            (z_ * vec.z_);
}
