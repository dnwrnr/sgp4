#ifndef SATELLITEORBIT_H_
#define SATELLITEORBIT_H_

#include "Tle.h"
#include "SGDP4.h"

class SatelliteOrbit {
public:
    SatelliteOrbit(void);
    virtual ~SatelliteOrbit(void);

    void SetTle(const Tle& tle);

    bool IsGeostationary();

    unsigned int GetOrbitNumber(const Julian& jul) const;

private:
    SGDP4 sgdp4_;
};

#endif

