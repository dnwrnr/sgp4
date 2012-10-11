#ifndef ECI_H_
#define ECI_H_

#include "CoordGeodetic.h"
#include "Vector.h"
#include "DateTime.h"

/**
 * A class store to store an Earth-centered inertial position.
 * This is only valid for the date specified.
 */
class Eci
{
public:

    /**
     * @param[in] dt the date to be used for this position
     * @param[in] latitude the latitude in degrees
     * @param[in] longitude the longitude in degrees
     * @param[in] altitude the altitude in kilometers
     */
    Eci(const DateTime& dt,
            const double latitude,
            const double longitude,
            const double altitude)
    {
        ToEci(dt, CoordGeodetic(latitude, longitude, altitude));
    }

    /**
     * @param[in] dt the date to be used for this position
     * @param[in] geo the position
     */
    Eci(const DateTime& dt, const CoordGeodetic& geo)
    {
        ToEci(dt, geo);
    }

    /**
     * @param[in] dt the date to be used for this position
     * @param[in] position
     */
    Eci(const DateTime &dt, const Vector &position)
        : m_dt(dt),
        m_position(position)
    {
    }

    /**
     * @param[in] dt the date to be used for this position
     * @param[in] position the position
     * @param[in] velocity the velocity
     */
    Eci(const DateTime &dt, const Vector &position, const Vector &velocity)
        : m_dt(dt),
        m_position(position),
        m_velocity(velocity)
    {
    }

    /**
     * Destructor
     */
    virtual ~Eci()
    {
    }

    /**
     * Equality operator
     * @param dt the date to compare
     * @returns true if the object matches
     */
    bool operator==(const DateTime& dt) const
    {
        return m_dt == dt;
    }

    /**
     * Inequality operator
     * @param dt the date to compare
     * @returns true if the object doesn't match
     */
    bool operator!=(const DateTime& dt) const
    {
        return m_dt != dt;
    }

    /**
     * Update this object with a new date and geodetic position
     * @param dt new date
     * @param geo new geodetic position
     */
    void Update(const DateTime& dt, const CoordGeodetic& geo)
    {
        ToEci(dt, geo);
    }

    /**
     * @returns the position
     */
    Vector GetPosition() const
    {
        return m_position;
    }

    /**
     * @returns the velocity
     */
    Vector GetVelocity() const
    {
        return m_velocity;
    }

    /**
     * @returns the date
     */
    DateTime GetDateTime() const
    {
        return m_dt;
    }

    /**
     * @returns the position in geodetic form
     */
    CoordGeodetic ToGeodetic() const;

private:
    void ToEci(const DateTime& dt, const CoordGeodetic& geo);

    DateTime m_dt;
    Vector m_position;
    Vector m_velocity;
};

#endif
