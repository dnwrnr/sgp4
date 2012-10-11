#ifndef OBSERVER_H_
#define OBSERVER_H_

#include "CoordGeodetic.h"
#include "Eci.h"

class DateTime;
class CoordTopographic;

class Observer
{
public:
    /**
     * Constructor
     * @param[in] latitude observers latitude in degrees
     * @param[in] longitude observers longitude in degrees
     * @param[in] altitude observers altitude in kilometers
     */
    Observer(const double latitude,
            const double longitude,
            const double altitude)
        : m_geo(latitude, longitude, altitude),
        m_eci(DateTime(), m_geo)
    {
    }

    /**
     * Constructor
     * @param[in] geo the observers position
     */
    Observer(const CoordGeodetic &geo)
        : m_geo(geo),
        m_eci(DateTime(), geo)
    {
    }

    /**
     * Destructor
     */
    virtual ~Observer()
    {
    }

    /**
     * Set the observers location
     * @param[in] geo the observers position
     */
    void SetLocation(const CoordGeodetic& geo)
    {
        m_geo = geo;
        m_eci.Update(m_eci.GetDateTime(), m_geo);
    }

    /**
     * Get the observers location
     * @returns the observers position
     */
    CoordGeodetic GetLocation() const
    {
        return m_geo;
    }

    /**
     * Get the look angle for the observers position to the object
     * @param[in] eci the object to find the look angle to
     * @returns the lookup angle
     */
    CoordTopographic GetLookAngle(const Eci &eci);

private:
    /**
     * @param[in] dt the date to update the observers position for
     */
    void Update(const DateTime &dt)
    {
        if (m_eci != dt)
        {
            m_eci.Update(dt, m_geo);
        }
    }

    /** the observers position */
    CoordGeodetic m_geo;
    /** the observers Eci for a particular time */
    Eci m_eci;
};

#endif

