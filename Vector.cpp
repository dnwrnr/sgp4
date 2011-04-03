#include "Vector.h"

double Vector::GetMagnitude() const {
    return sqrt(x_ * x_ + y_ * y_ + z_ * z_);
}
