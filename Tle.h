#ifndef TLE_H_
#define TLE_H_

#include "Globals.h"
#include "Julian.h"

#include <iostream>

class Tle {
public:
    Tle(const std::string& line_one, const std::string& line_two);
    Tle(const std::string& name, const std::string& line_one, const std::string& line_two);
    Tle(const Tle& tle);
    virtual ~Tle();

    /*
     * get raw strings
     */
    std::string GetName() const {
        return name_;
    }

    std::string GetLine1() const {
        return line_one_;
    }

    std::string GetLine2() const {
        return line_two_;
    }

    /*
     * get tle values
     */
    unsigned int NoradNumber() const {
        return norad_number_;
    }

    std::string InternationlDesignator() const {
        return international_designator_;
    }

    Julian Epoch()const {
        return epoch_;
    }

    double MeanMotionDot()const {
        return mean_motion_dot_;
    }

    double MeanMotionDot2()const {
        return mean_motion_dot2_;
    }

    double BStar()const {
        return bstar_;
    }

    double Inclination(const bool in_degrees)const {
        if (in_degrees)
            return inclination_;
        else
            return DegreesToRadians(inclination_);
    }

    double RightAscendingNode(const bool in_degrees)const {
        if (in_degrees)
            return right_ascending_node_;
        else
            return DegreesToRadians(right_ascending_node_);
    }

    double Eccentricity()const {
        return eccentricity_;
    }

    double ArgumentPerigee(const bool in_degrees)const {
        if (in_degrees)
            return argument_perigee_;
        else
            return DegreesToRadians(argument_perigee_);
    }

    double MeanAnomaly(const bool in_degrees)const {
        if (in_degrees)
            return mean_anomaly_;
        else
            return DegreesToRadians(mean_anomaly_);
    }

    double MeanMotion()const {
        return mean_motion_;
    }

    unsigned int OrbitNumber()const {
        return orbit_number_;
    }

    /*
     * helper / validation methods
     */
    static unsigned int GetLineLength() {
        return TLE_LEN_LINE_DATA;
    }
    static bool IsValidPair(const std::string& line1, const std::string& line2);
    static bool IsValidLine(const std::string& str, unsigned char line_number);

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
    static bool ValidateLine(const std::string& line, const std::string& pattern);
    /*
     * compute checksum
     */
    static int CheckSum(const std::string& str);

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
     * name line
     */
    static const unsigned int TLE_LEN_LINE_DATA = 69;
    static const unsigned int TLE_LEN_LINE_NAME = 22;

    /*
     * line 1
     */
    static const unsigned int TLE1_COL_NORADNUM = 2;
    static const unsigned int TLE1_LEN_NORADNUM = 5;
    static const unsigned int TLE1_COL_INTLDESC_A = 9;
    static const unsigned int TLE1_LEN_INTLDESC_A = 2;
    static const unsigned int TLE1_COL_INTLDESC_B = 11;
    static const unsigned int TLE1_LEN_INTLDESC_B = 3;
    static const unsigned int TLE1_COL_INTLDESC_C = 14;
    static const unsigned int TLE1_LEN_INTLDESC_C = 3;
    static const unsigned int TLE1_COL_EPOCH_A = 18;
    static const unsigned int TLE1_LEN_EPOCH_A = 2;
    static const unsigned int TLE1_COL_EPOCH_B = 20;
    static const unsigned int TLE1_LEN_EPOCH_B = 12;
    static const unsigned int TLE1_COL_MEANMOTIONDT = 33;
    static const unsigned int TLE1_LEN_MEANMOTIONDT = 10;
    static const unsigned int TLE1_COL_MEANMOTIONDT2 = 44;
    static const unsigned int TLE1_LEN_MEANMOTIONDT2 = 8;
    static const unsigned int TLE1_COL_BSTAR = 53;
    static const unsigned int TLE1_LEN_BSTAR = 8;
    static const unsigned int TLE1_COL_EPHEMTYPE = 62;
    static const unsigned int TLE1_LEN_EPHEMTYPE = 1;
    static const unsigned int TLE1_COL_ELNUM = 64;
    static const unsigned int TLE1_LEN_ELNUM = 4;

    /*
     * line 2
     */
    static const unsigned int TLE2_COL_NORADNUM = 2;
    static const unsigned int TLE2_LEN_NORADNUM = 5;
    static const unsigned int TLE2_COL_INCLINATION = 8;
    static const unsigned int TLE2_LEN_INCLINATION = 8;
    static const unsigned int TLE2_COL_RAASCENDNODE = 17;
    static const unsigned int TLE2_LEN_RAASCENDNODE = 8;
    static const unsigned int TLE2_COL_ECCENTRICITY = 26;
    static const unsigned int TLE2_LEN_ECCENTRICITY = 7;
    static const unsigned int TLE2_COL_ARGPERIGEE = 34;
    static const unsigned int TLE2_LEN_ARGPERIGEE = 8;
    static const unsigned int TLE2_COL_MEANANOMALY = 43;
    static const unsigned int TLE2_LEN_MEANANOMALY = 8;
    static const unsigned int TLE2_COL_MEANMOTION = 52;
    static const unsigned int TLE2_LEN_MEANMOTION = 11;
    static const unsigned int TLE2_COL_REVATEPOCH = 63;
    static const unsigned int TLE2_LEN_REVATEPOCH = 5;
};

#endif


