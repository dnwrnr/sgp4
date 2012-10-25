#ifndef SOLARPOSITION_H_
#define SOLARPOSITION_H_

#include "DateTime.h"
#include "Eci.h"

/**
 * @brief Find the position of the sun
 */
class SolarPosition
{
public:
    SolarPosition()
    {
    }

    virtual ~SolarPosition()
    {
    }

    Eci FindPosition(const DateTime& dt);

private:
    double Delta_ET(double year) const;
};

#endif
