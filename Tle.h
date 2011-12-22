#ifndef TLE_H_
#define TLE_H_

#include "Globals.h"
#include "Util.h"
#include "Julian.h"
#include "TleException.h"

#include <iostream>

class Tle
{
public:
    Tle(const std::string& line_one, const std::string& line_two)
    : line_one_(line_one), line_two_(line_two)
    {
        Initialize();
    }

    Tle(const std::string& name, const std::string& line_one,
        const std::string& line_two)
    : name_(name),line_one_(line_one), line_two_(line_two)
    {
        Initialize();
    }

    Tle(const Tle& tle);

    virtual ~Tle()
    {
    }

    /*
     * get raw strings
     */
    std::string GetName() const
    {
        return name_;
    }

    std::string GetLine1() const
    {
        return line_one_;
    }

    std::string GetLine2() const
    {
        return line_two_;
    }

    /*
     * get tle values
     */
    unsigned int NoradNumber() const
    {
        return norad_number_;
    }

    std::string InternationlDesignator() const
    {
        return international_designator_;
    }

    Julian Epoch() const
    {
        return epoch_;
    }

    double MeanMotionDot() const
    {
        return mean_motion_dot_;
    }

    double MeanMotionDot2() const
    {
        return mean_motion_dot2_;
    }

    double BStar() const
    {
        return bstar_;
    }

    double Inclination(bool in_degrees) const
    {
        if (in_degrees)
        {
            return inclination_;
        }
        else
        {
            return Util::DegreesToRadians(inclination_);
        }
    }

    double RightAscendingNode(const bool in_degrees) const
    {
        if (in_degrees)
        {
            return right_ascending_node_;
        }
        else
        {
            return Util::DegreesToRadians(right_ascending_node_);
        }
    }

    double Eccentricity() const
    {
        return eccentricity_;
    }

    double ArgumentPerigee(const bool in_degrees) const
    {
        if (in_degrees)
        {
            return argument_perigee_;
        }
        else
        {
            return Util::DegreesToRadians(argument_perigee_);
        }
    }

    double MeanAnomaly(const bool in_degrees) const
    {
        if (in_degrees)
        {
            return mean_anomaly_;
        }
        else
        {
            return Util::DegreesToRadians(mean_anomaly_);
        }
    }

    double MeanMotion() const
    {
        return mean_motion_;
    }

    unsigned int OrbitNumber() const
    {
        return orbit_number_;
    }

    /*
     * helper / validation methods
     */
    static unsigned int GetLineLength()
    {
        return TLE_LEN_LINE_DATA;
    }

    static void IsValidPair(const std::string& line1, const std::string& line2);
    static void IsValidLine(const std::string& str, int line_number);

private:
    /*
     * initialize from raw tle strings
     */
    void Initialize();
    /*
     * format a raw string into an exponent string
     */
    static std::string ExpToDecimal(const std::string&);
    /*
     * validate a line against a pattern
     */
    static void ValidateLine(const std::string& line,
            const std::string& pattern);
    /*
     * compute checksum
     */
    static int CheckSum(const std::string& str);

    static std::string ExtractNoradNumber(const std::string& str,
            int line_number);

    static bool IsValidLineLength(const std::string& str);

private:
    /*
     * raw tle data
     */
    std::string name_;
    std::string line_one_;
    std::string line_two_;

    /*
     * extracted values all in native units
     */
    unsigned int norad_number_;
    std::string international_designator_;
    Julian epoch_;
    double mean_motion_dot_;
    double mean_motion_dot2_;
    double bstar_;
    double inclination_;
    double right_ascending_node_;
    double eccentricity_;
    double argument_perigee_;
    double mean_anomaly_;
    double mean_motion_;
    unsigned int orbit_number_;

    /*
     * line lengths
     */
    static const unsigned int TLE_LEN_LINE_DATA = 69;
    static const unsigned int TLE_LEN_LINE_NAME = 22;

    };

#endif

