#ifndef ECI_H_
#define ECI_H_

#include "Coord.h"
#include "Vector.h"
#include "Julian.h"

class Eci {
public:

    Eci() {
    };
    Eci(const Julian &date, const CoordGeodetic &geo);
    Eci(const Julian &date, const Vector &position, const Vector &velocity);
    virtual ~Eci(void);

    Vector GetPosition() const {
        return position_;
    }

    Vector GetVelocity() const {
        return velocity_;
    }

    Julian GetDate() const {
        return date_;
    }

private:
    Vector position_;
    Vector velocity_;
    Julian date_;
};

#endif

