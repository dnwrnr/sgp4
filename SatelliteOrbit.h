#ifndef SATELLITEORBIT_H_
#define SATELLITEORBIT_H_

#include "Tle.h"
#include "SGP4.h"

class SatelliteOrbit {
public:
    SatelliteOrbit(void);
    virtual ~SatelliteOrbit(void);

    void SetTle(const Tle& tle);

    bool IsGeostationary();

    unsigned int GetOrbitNumber(const Julian& jul) const;

    static double Footprint(const double& altitude);

private:
    SGP4 sgp4_;
};

#endif

