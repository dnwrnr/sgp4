#ifndef VECTOR_H_
#define VECTOR_H_

#include <cmath>

struct Vector {
public:

    Vector(void)
    : x(0.0), y(0.0), z(0.0), w(0.0) {
    }

    Vector(double x_in, double y_in, double z_in)
    : x(x_in), y(y_in), z(z_in), w(0.0) {
    }

    Vector(double x_in, double y_in, double z_in, double w_in)
    : x(x_in), y(y_in), z(z_in), w(w_in) {
    }

    virtual ~Vector() {
    };

    double GetMagnitude() const;
    Vector Subtract(const Vector& vec) const;
    double Dot(const Vector& vec) const;

    double x;
    double y;
    double z;
    double w;
};

#endif
