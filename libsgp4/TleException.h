#ifndef TLEEXCEPTION_H_
#define TLEEXCEPTION_H_

#include <exception>
#include <iostream>

class TleException : public std::exception
{
public:

    TleException(const char* message)
        : message_(message)
    {
    }

    virtual ~TleException(void) throw ()
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
