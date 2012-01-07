#ifndef DECAYEDEXCEPTION_H_
#define DECAYEDEXCEPTION_H_

#include <exception>
#include <iostream>
#include "Julian.h"
#include "Vector.h"

class DecayedException : public std::exception
{
public:
    DecayedException(const Julian& dt, const Vector& pos, const Vector& vel)
        : _dt(dt), _pos(pos), _vel(vel)
    {
    }

    virtual ~DecayedException(void) throw ()
    {
    }

    virtual const char* what() const throw ()
    {
        return "Error: Satellite decayed";
    }

    Julian GetDate() const
    {
        return _dt;
    }

    Vector GetPosition() const
    {
        return _pos;
    }

    Vector GetVelocity() const
    {
        return _vel;
    }

private:
    Julian _dt;
    Vector _pos;
    Vector _vel;
};

#endif
