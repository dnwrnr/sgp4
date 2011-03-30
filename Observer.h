#ifndef OBSERVER_H_
#define OBSERVER_H_

#include "Coord.h"
#include "Julian.h"
#include "Eci.h"

class Observer {
public:
    Observer(const double latitude, const double longitude, const double altitude);
    Observer(const CoordGeodetic &geo);
    virtual ~Observer(void);

    void SetLocation(const CoordGeodetic& geo) {
        geo_ = geo;
        observers_eci_ = Eci(observers_eci_.GetDate(), geo_);
    }

    CoordGeodetic GetLocation() const {
        return geo_;
    }

    Eci GetEciPosition(const Julian &date) const {
        return Eci(date, geo_);
    }

    CoordTopographic GetLookAngle(const Eci &eci);

private:
    void UpdateObserversEci(const Julian &date);

    /*
     * the observers eci for a particular time
     */
    Eci observers_eci_;

    /*
     * the observers position
     */
    CoordGeodetic geo_;
};

#endif

