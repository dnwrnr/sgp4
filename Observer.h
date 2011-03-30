#ifndef OBSERVER_H_
#define OBSERVER_H_

#include "Coord.h"

class Observer {
public:
    Observer(const double latitude, const double longitude, const double altitude);
    Observer(const CoordGeodetic &geo);
    virtual ~Observer(void);

    void SetLocation(const CoordGeodetic& geo) {
        geo_ = geo;
    }

    CoordGeodetic GetLocation() const {
        return geo_;
    }

private:
    /*
     * the observers position
     */
    CoordGeodetic geo_;
};

#endif

