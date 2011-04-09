#ifndef VECTOR_H_
#define VECTOR_H_

#include <cmath>

class Vector {
public:

    Vector(void)
    : x_(0.0), y_(0.0), z_(0.0), w_(0.0) {
    }

    Vector(double x, double y, double z)
    : x_(x), y_(y), z_(z), w_(0.0) {
    }

    Vector(double x, double y, double z, double w)
    : x_(x), y_(y), z_(z), w_(w) {
    }

    virtual ~Vector() {
    };

    void SetX(const double& x) {
        x_ = x;
    }

    void SetY(const double& y) {
        y_ = y;
    }

    void SetZ(const double& z) {
        z_ = z;
    }

    void SetW(const double& w) {
        w_ = w;
    }

    double GetX() const {
        return x_;
    }

    double GetY() const {
        return y_;
    }

    double GetZ() const {
        return z_;
    }

    double GetW() const {
        return w_;
    }

    double GetMagnitude() const;
    Vector Subtract(const Vector& vec) const;
    double Dot(const Vector& vec) const;

protected:
    double x_;
    double y_;
    double z_;
    double w_;
};

#endif