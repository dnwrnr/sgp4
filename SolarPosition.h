#ifndef SOLARPOSITION_H_
#define SOLARPOSITION_H_

#include "Julian.h"
#include "Eci.h"

class SolarPosition {
public:
    SolarPosition(void);
    virtual ~SolarPosition(void);

    void FindPosition(const Julian& j, Eci& eci);

private:
    double Modulus(double arg1, double arg2) const;
    double Delta_ET(double year) const;
};

#endif
