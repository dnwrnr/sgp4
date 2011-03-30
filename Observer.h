#ifndef OBSERVER_H_
#define OBSERVER_H_

#include "Coord.h"

class Observer {
public:
    Observer(void);
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

