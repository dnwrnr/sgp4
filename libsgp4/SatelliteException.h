#ifndef SATELLITEEXCEPTION_H_
#define SATELLITEEXCEPTION_H_

#include <exception>

/**
 * @brief The exception that the SGP4 class throws upon an error.
 */
class SatelliteException : public std::exception
{
public:
    SatelliteException(const char* message)
        : message_(message)
    {
    }

    virtual ~SatelliteException(void) throw ()
    {
    }

    virtual const char* what() const throw ()
    {
        return message_.c_str();
    }

private:
    std::string message_;
};

#endif
