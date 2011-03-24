#ifndef VECTOR_H_
#define VECTOR_H_

class Vector {
public:

    Vector(double x = 0.0, double y = 0.0, double z = 0.0, double w = 0.0)
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

protected:
    double x_;
    double y_;
    double z_;
    double w_;
};

#endif