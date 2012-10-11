#ifndef TLEEXCEPTION_H_
#define TLEEXCEPTION_H_

#include <exception>

class TleException : public std::exception
{
public:
    /**
     * Constructor
     * @param message Exception message
     */
    TleException(const char* message)
        : m_message(message)
    {
    }

    /**
     * Destructor
     */
    virtual ~TleException(void) throw ()
    {
    }

    /**
     * Get the exception message
     * @returns the exception message
     */
    virtual const char* what() const throw ()
    {
        return m_message.c_str();
    }

private:
    /** the exception message */
    std::string m_message;
};

#endif
