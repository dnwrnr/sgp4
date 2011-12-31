#ifndef SOLARPOSITION_H_
#define SOLARPOSITION_H_

#include "Julian.h"
#include "Eci.h"

class SolarPosition
{
public:
    SolarPosition()
    {
    }

    virtual ~SolarPosition()
    {
    }

    Eci FindPosition(const Julian& j);

private:
    double Delta_ET(double year) const;
};

#endif
