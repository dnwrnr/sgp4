#ifndef ECI_H_
#define ECI_H_

#include "CoordGeodetic.h"
#include "Vector.h"
#include "Julian.h"
#include "Globals.h"

class Eci
{
public:

    /*
     * in degrees
     */
    Eci(const Julian& date, double latitude, double longitude, double altitude)
    {
        ToEci(date, CoordGeodetic(latitude, longitude, altitude));
    }

    Eci(const Julian& date, const CoordGeodetic& g)
    {
        ToEci(date, g);
    }

    Eci(const Julian &date, const Vector &position)
        : date_(date), position_(position)
    {
    }

    Eci(const Julian &date, const Vector &position, const Vector &velocity)
        : date_(date), position_(position), velocity_(velocity)
    {
    }

    virtual ~Eci()
    {
    }

    Vector GetPosition() const
    {
        return position_;
    }

    Vector GetVelocity() const
    {
        return velocity_;
    }

    Julian GetDate() const
    {
        return date_;
    }

    CoordGeodetic ToGeodetic() const;

protected:
    void ToEci(const Julian& date, double latitude, double longitude,
            double altitude);
    void ToEci(const Julian& date, const CoordGeodetic& g);

private:
    Julian date_;
    Vector position_;
    Vector velocity_;
};

#endif

