#ifndef OBSERVER_H_
#define OBSERVER_H_

#include "CoordGeodetic.h"
#include "CoordTopographic.h"
#include "Julian.h"
#include "Eci.h"

class Observer
{
public:
    Observer(double latitude, double longitude, double altitude)
        : geo_(latitude, longitude, altitude),
        observers_eci_(Julian(), geo_)
    {
    }

    Observer(const CoordGeodetic &g)
        : geo_(g), observers_eci_(Julian(), geo_)
    {
    }

    virtual ~Observer()
    {
    }

    void SetLocation(const CoordGeodetic& g)
    {
        geo_ = g;
        observers_eci_ = Eci(observers_eci_.GetDate(), geo_);
    }

    CoordGeodetic GetLocation() const
    {
        return geo_;
    }

    Eci GetEciPosition(const Julian &date) const
    {
        return Eci(date, geo_);
    }

    CoordTopographic GetLookAngle(const Eci &eci);

private:
    void UpdateObserversEci(const Julian &date)
    {
        /*
         * if date has changed, update for new date
         */
        if (observers_eci_.GetDate() != date)
        {
            observers_eci_ = Eci(date, geo_);
        }
    }

    /*
     * the observers position
     */
    CoordGeodetic geo_;
    /*
     * the observers eci for a particular time
     */
    Eci observers_eci_;
};

#endif

