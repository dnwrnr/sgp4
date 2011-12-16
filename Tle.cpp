#include "Tle.h"

#include <cstdlib>

namespace
{
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

}

Tle::Tle(const Tle& tle)
{
    name_ = tle.name_;
    line_one_ = tle.line_one_;
    line_two_ = tle.line_two_;

    norad_number_ = tle.norad_number_;
    international_designator_ = tle.international_designator_;
    epoch_ = tle.epoch_;
    mean_motion_dot_ = tle.mean_motion_dot_;
    mean_motion_dot2_ = tle.mean_motion_dot2_;
    bstar_ = tle.bstar_;
    inclination_ = tle.inclination_;
    right_ascending_node_ = tle.right_ascending_node_;
    eccentricity_ = tle.eccentricity_;
    argument_perigee_ = tle.argument_perigee_;
    mean_anomaly_ = tle.mean_anomaly_;
    mean_motion_ = tle.mean_motion_;
    orbit_number_ = tle.orbit_number_;
}

/*
 * convert a tle raw string into an exponent string
 */
std::string Tle::ExpToDecimal(const std::string& str)
{
    static const int START_SIGN = 0;
    static const int LENGTH_SIGN = 1;
    static const int START_MANTISSA = 1;
    static const int LENGTH_MANTISSA = 5;
    static const int START_EXP = 6;
    static const int LENGTH_EXP = 2;

    if ((LENGTH_SIGN + LENGTH_MANTISSA + LENGTH_EXP) != str.length())
    {
        throw TleException("Invalid string length for exponential conversion.");
    }

    std::string value = str.substr(START_SIGN, LENGTH_SIGN);
    value += ".";
    value += str.substr(START_MANTISSA, LENGTH_MANTISSA);
    value += "e";
    value += str.substr(START_EXP, LENGTH_EXP);

    return value;
}

/*
 * extract all variables
 */
void Tle::Initialize()
{
    std::string temp;

    /*
     * trim whitespace
     */
    TrimLeft(name_);
    TrimRight(name_);
    TrimLeft(line_one_);
    TrimRight(line_one_);
    TrimLeft(line_two_);
    TrimRight(line_two_);

    /*
     * check the two lines are valid
     */
    IsValidPair(line_one_, line_two_);

    /*
     * line 1
     */

    temp = ExtractNoradNumber(line_one_, 1);
    norad_number_ = atoi(temp.c_str());
    /*
     * if blank use norad number for name
     */
    if (name_.empty())
    {
        name_ = temp;
    }

    international_designator_ = line_one_.substr(TLE1_COL_INTLDESC_A,
            TLE1_LEN_INTLDESC_A + TLE1_LEN_INTLDESC_B + TLE1_LEN_INTLDESC_C);

    int year = atoi(line_one_.substr(TLE1_COL_EPOCH_A,
                TLE1_LEN_EPOCH_A).c_str());
    double day = atof(line_one_.substr(TLE1_COL_EPOCH_B,
                TLE1_LEN_EPOCH_B).c_str());
    /*
     * generate julian date for epoch
     */
    if (year < 57)
        year += 2000;
    else
        year += 1900;
    epoch_ = Julian(year, day);

    if (line_one_[TLE1_COL_MEANMOTIONDT] == '-') {
        temp = "-0";
    } else
        temp = "0";
    temp += line_one_.substr(TLE1_COL_MEANMOTIONDT + 1, TLE1_LEN_MEANMOTIONDT);
    mean_motion_dot_ = atof(temp.c_str());

    temp = ExpToDecimal(line_one_.substr(TLE1_COL_MEANMOTIONDT2,
                TLE1_LEN_MEANMOTIONDT2));
    mean_motion_dot2_ = atof(temp.c_str());

    temp = ExpToDecimal(line_one_.substr(TLE1_COL_BSTAR,
                TLE1_LEN_BSTAR).c_str());
    bstar_ = atof(temp.c_str());

    /*
     * line 2
     */

    temp = line_two_.substr(TLE2_COL_INCLINATION, TLE2_LEN_INCLINATION);
    TrimLeft(temp);
    inclination_ = atof(temp.c_str());

    temp = line_two_.substr(TLE2_COL_RAASCENDNODE, TLE2_LEN_RAASCENDNODE);
    TrimLeft(temp);
    right_ascending_node_ = atof(temp.c_str());

    temp = "0.";
    temp += line_two_.substr(TLE2_COL_ECCENTRICITY, TLE2_LEN_ECCENTRICITY);
    eccentricity_ = atof(temp.c_str());

    temp = line_two_.substr(TLE2_COL_ARGPERIGEE, TLE2_LEN_ARGPERIGEE);
    TrimLeft(temp);
    argument_perigee_ = atof(temp.c_str());

    temp = line_two_.substr(TLE2_COL_MEANANOMALY, TLE2_LEN_MEANANOMALY);
    TrimLeft(temp);
    mean_anomaly_ = atof(temp.c_str());

    temp = line_two_.substr(TLE2_COL_MEANMOTION, TLE2_LEN_MEANMOTION);
    TrimLeft(temp);
    mean_motion_ = atof(temp.c_str());

    temp = line_two_.substr(TLE2_COL_REVATEPOCH, TLE2_LEN_REVATEPOCH);
    TrimLeft(temp);
    orbit_number_ = atoi(temp.c_str());
}

/*
 * check the two lines have matching norad numbers
 * and that the lines themselves are equal
 */
void Tle::IsValidPair(const std::string& line1, const std::string& line2)
{
    /*
     * validate each line
     */
    IsValidLine(line1, 1);
    IsValidLine(line2, 2);

    /*
     * extract norad numbers
     */
    std::string norad_1 = ExtractNoradNumber(line1, 1);
    std::string norad_2 = ExtractNoradNumber(line2, 2);

    /*
     * make sure they match
     */
    if (norad_1.compare(norad_2) != 0)
    {
        throw TleException("Norad numbers do not match.");
    }
}

/*
 * validate a line
 */
void Tle::IsValidLine(const std::string& str, int line_number)
{
    /*
     * validation patterns
     */
    static const std::string line1_pattern = "1 NNNNNC NNNNNXXX NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN";
    static const std::string line2_pattern = "2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN";

    /*
     * validate line against the pattern
     */
    if (1 == line_number)
    {
        ValidateLine(str, line1_pattern);
    }
    else if (2 == line_number)
    {
        ValidateLine(str, line2_pattern);
    }
    else
    {
        throw TleException("Invalid line number to check.");
    }

    /*
     * last char in string is modulo 10 checksum
     * edited out as checksum isnt consistent
     *
     * int chk = CheckSum(str);
     * if (chk != (str[TLE_LEN_LINE_DATA - 1] - '0'))
     *   return false;
     */
}

bool Tle::IsValidLineLength(const std::string& str)
{
    return str.length() == GetLineLength() ? true : false;
}

/*
 * validate line given a pattern
 */
void Tle::ValidateLine(const std::string& line, const std::string& pattern)
{
    /*
     * check length of line
     */
    if (!IsValidLineLength(line))
    {
        throw TleException("Invalid line length.");
    }

    std::string::const_iterator pattern_itr = pattern.begin();
    std::string::const_iterator line_itr = line.begin();

    while (pattern_itr != pattern.end())
    {
        if (isdigit(*pattern_itr) || *pattern_itr == ' ' ||
                *pattern_itr == '.')
        {
            /*
             * should match exactly
             */
            if (*pattern_itr != *line_itr)
            {
                throw TleException("Invalid character.");
            }

        }
        else if (*pattern_itr == 'N')
        {
            /*
             * 'N' = number or ' '
             */
            if (!isdigit(*line_itr) && *line_itr != ' ')
            {
                throw TleException("Invalid character.");
            }

        }
        else if (*pattern_itr == '+')
        {
            /*
             * '+' = '+' or '-' or ' ' or '0'
             */
            if (*line_itr != '+' && *line_itr != '-' &&
                    *line_itr != ' ' && *line_itr != '0')
            {
                throw TleException("Invalid character.");
            }

        }
        else if (*pattern_itr == '-')
        {
            /*
             * '-' = '+' or '-'
             */
            if (*line_itr != '+' && *line_itr != '-')
            {
                throw TleException("Invalid character.");
            }

        }
        else if (*pattern_itr == 'C')
        {
            /*
             * 'C' = 'U' or 'S'
             */
            if (*line_itr != 'U' && *line_itr != 'S')
            {
                throw TleException("Invalid character.");
            }

        }
        else if (*pattern_itr == 'X')
        {
            /*
             * 'X' = A-Z or ' '
             */
            if (!(*line_itr >= 'A' || *line_itr <= 'Z') && *line_itr != ' ')
            {
                throw TleException("Invalid character.");
            }
        }

        pattern_itr++;
        line_itr++;
    }
}

/*
 * compute checksum
 */
int Tle::CheckSum(const std::string & str)
{
    size_t len = str.size() - 1;
    int xsum = 0;

    for (size_t i = 0; i < len; i++)
    {
        char ch = str[i];

        if (isdigit(ch))
        {
            xsum += (ch - '0');
        }
        else if (ch == '-')
        {
            xsum++;
        }
    }

    return (xsum % 10);
}

std::string Tle::ExtractNoradNumber(const std::string& str, int line_number)
{
    std::string norad_number;

    /*
     * check length
     */
    if (!IsValidLineLength(str))
    {
        throw TleException("Invalid line length.");
    }

    /*
     * extract string
     */
    if (1 == line_number)
    {
        norad_number = str.substr(TLE1_COL_NORADNUM, TLE1_LEN_NORADNUM);
    }
    else if (2 == line_number)
    {
        norad_number = str.substr(TLE2_COL_NORADNUM, TLE2_LEN_NORADNUM);
    }
    else
    {
        throw TleException("Invalid line number to check.");
    }

    return norad_number;
}
