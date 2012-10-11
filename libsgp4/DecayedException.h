#ifndef DECAYEDEXCEPTION_H_
#define DECAYEDEXCEPTION_H_

#include <exception>

#include "DateTime.h"
#include "Vector.h"

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
    DateTime GetDateTime() const
    {
        return _dt;
    }

    /**
     * @returns the position
     */
    Vector GetPosition() const
    {
        return _pos;
    }

    /**
     * @returns the velocity
     */
    Vector GetVelocity() const
    {
        return _vel;
    }

private:
    DateTime _dt;
    Vector _pos;
    Vector _vel;
};

#endif
