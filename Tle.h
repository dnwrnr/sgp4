#ifndef TLE_H_
#define TLE_H_

#include "Globals.h"
#include "Julian.h"

#include <iostream>

class Tle {
public:
    Tle(const std::string& name, const std::string& line_one, const std::string& line_two);
    Tle(const Tle& tle);
    virtual ~Tle();

    enum TleLine {
        LINE_ZERO,
        LINE_ONE,
        LINE_TWO
    };

    enum TleField {
        FLD_FIRST,
        FLD_NORADNUM = FLD_FIRST,
        FLD_INTLDESC,
        FLD_SET, // Tle set number
        FLD_EPOCHYEAR, // Epoch: Last two digits of year
        FLD_EPOCHDAY, // Epoch: Fractional Julian Day of year
        FLD_ORBITNUM, // Orbit at epoch
        FLD_I, // Inclination
        FLD_RAAN, // R.A. ascending node
        FLD_E, // Eccentricity
        FLD_ARGPER, // Argument of perigee
        FLD_M, // Mean anomaly
        FLD_MMOTION, // Mean motion
        FLD_MMOTIONDT, // First time derivative of mean motion
        FLD_MMOTIONDT2, // Second time derivative of mean motion
        FLD_BSTAR, // BSTAR Drag
        FLD_LAST // MUST be last
    };

    enum TleUnits {
        U_FIRST,
        U_RAD = U_FIRST, // radians
        U_DEG, // degrees
        U_NATIVE, // Tle format native units (no conversion)
        U_LAST // MUST be last
    };

    /*
     * initialize from raw tle strings
     */
    void Initialize();

    /*
     * get field represented as a double
     */
    double GetField(TleField fld, TleUnits unit = U_NATIVE) const;
    /*
     * get field represented as a string
     */
    std::string GetFieldString(TleField fld, bool append_units) const;

    /*
     * get epoch of tle
     */
    Julian GetEpoch() const {
        return date_;
    }

    /*
     * compute checksum
     */
    static int CheckSum(const std::string& str);
    static bool IsValidLine(std::string& str, TleLine line);
    static void TrimLeft(std::string& str);
    static void TrimRight(std::string& str);

    std::string GetName() const {
        return name_;
    }

    std::string GetLine1() const {
        return line_one_;
    }

    std::string GetLine2() const {
        return line_two_;
    }

    double GetNoradNumber() const {
        return GetField(Tle::FLD_NORADNUM);
    }

private:
    /*
     * get units for a field
     */
    std::string GetUnits(TleField field) const;
    /*
     * get raw value for a field
     */
    double GetFieldNumeric(TleField) const;
    /*
     * convert a raw string into an exponent string
     */
    static std::string ExpToDecimal(const std::string&);
    /*
     * validate a line
     */
    static bool ValidateLine(const std::string& pattern, const std::string& line);

    /*
     * raw tle data
     */
    std::string name_;
    std::string line_one_;
    std::string line_two_;

    /*
     * epoch of tle
     */
    Julian date_;

    // array of tle values in native string and double form
    std::pair<std::string, double> fields_[FLD_LAST];

    /*
     * name line
     */
    static const unsigned int TLE_LEN_LINE_DATA = 69;
    static const unsigned int TLE_LEN_LINE_NAME = 22;

    /*
     * line 1
     */
    static const unsigned int TLE1_COL_SATNUM = 2;
    static const unsigned int TLE1_LEN_SATNUM = 5;
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
    static const unsigned int TLE2_COL_SATNUM = 2;
    static const unsigned int TLE2_LEN_SATNUM = 5;
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


