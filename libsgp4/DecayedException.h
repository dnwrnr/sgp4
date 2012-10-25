#ifndef DECAYEDEXCEPTION_H_
#define DECAYEDEXCEPTION_H_

#include "DateTime.h"
#include "Vector.h"

#include <exception>

/**
 * @brief The exception that the SGP4 class throws when a satellite decays.
 */
class DecayedException : public std::exception
{
public:
    /**
     * Constructor
     * @param[in] dt time of the event
     * @param[in] pos position of the satellite at dt
     * @param[in] vel velocity of the satellite at dt
     */
    DecayedException(const DateTime& dt, const Vector& pos, const Vector& vel)
        : _dt(dt), _pos(pos), _vel(vel)
    {
    }

    /**
     * Destructor
     */
    virtual ~DecayedException(void) throw ()
    {
    }

    /**
     * @returns the error string
     */
    virtual const char* what() const throw ()
    {
        return "Error: Satellite decayed";
    }

    /**
     * @returns the date
     */
    DateTime Decayed() const
    {
        return _dt;
    }

    /**
     * @returns the position
     */
    Vector Position() const
    {
        return _pos;
    }

    /**
     * @returns the velocity
     */
    Vector Velocity() const
    {
        return _vel;
    }

private:
    DateTime _dt;
    Vector _pos;
    Vector _vel;
};

#endif
