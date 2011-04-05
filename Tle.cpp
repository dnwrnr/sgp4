#include "Tle.h"

#include <cassert>
#include <cstdlib>

Tle::Tle(const std::string& name, const std::string& line_one, const std::string& line_two) {
    name_ = name;
    line_one_ = line_one;
    line_two_ = line_two;

    Initialize();
}

Tle::Tle(const Tle& tle) {
    name_ = tle.name_;
    line_one_ = tle.line_one_;
    line_two_ = tle.line_two_;
    date_ = tle.date_;

    for (int fld = FLD_FIRST; fld < FLD_LAST; fld++) {
        fields_[fld] = tle.fields_[fld];
    }
}

Tle::~Tle() {
}

double Tle::GetField(TleField fld, TleUnits units) const {
    assert((FLD_FIRST <= fld) && (fld < FLD_LAST));
    assert((U_FIRST <= units) && (units < U_LAST));

    double value = fields_[fld].second;
    double converted = 0.0;

    switch (fld) {
        case FLD_I:
        case FLD_RAAN:
        case FLD_ARGPER:
        case FLD_M:
        {
            /*
             * native format is degrees
             */
            if (units == U_RAD)
                converted = Globals::Deg2Rad(value);
            break;
        }
        default:
        {
            converted = value;
            break;
        }
    }

    return converted;
}

/*
 * get field represented as a string
 */
std::string Tle::GetFieldString(TleField fld, bool append_units) const {
    assert((FLD_FIRST <= fld) && (fld < FLD_LAST));

    std::string str = fields_[fld].first;

    if (append_units)
        str += GetUnits(fld);

    return str;
}

/*
 * get the units for a field
 */
std::string Tle::GetUnits(TleField field) const {
    static const std::string strDegrees = " degrees";
    static const std::string strRevsPerDay = " revs / day";
    static const std::string strNull = "";

    switch (field) {
        case FLD_I:
        case FLD_RAAN:
        case FLD_ARGPER:
        case FLD_M:
            return strDegrees;

        case FLD_MMOTION:
            return strRevsPerDay;

        default:
            return strNull;
    }
}

/*
 * convert a raw string into an exponent string
 */
std::string Tle::ExpToDecimal(const std::string& str) {
    static const int kColumnSign = 0;
    static const int kLengthSign = 1;
    static const int kColumnMantissa = 1;
    static const int kLengthMantissa = 5;
    static const int kColumnExp = 6;
    static const int kLengthExp = 2;

    assert(8 == str.length());

    std::string value;
    value = str.substr(kColumnSign, kLengthSign);
    value += ".";
    value += str.substr(kColumnMantissa, kLengthMantissa);
    value += "e";
    value += str.substr(kColumnExp, kLengthExp);

    return value;
}

/*
 * extract all variables
 */
void Tle::Initialize() {
    TrimLeft(name_);
    TrimRight(name_);
    TrimLeft(line_one_);
    TrimRight(line_one_);
    TrimLeft(line_two_);
    TrimRight(line_two_);

    assert(!name_.empty());
    assert(!line_one_.empty());
    assert(!line_two_.empty());

    assert(IsValidLine(line_one_, Tle::LINE_ONE));
    assert(IsValidLine(line_two_, Tle::LINE_TWO));

    /*
     * line 1
     */

    fields_[FLD_NORADNUM].first = line_one_.substr(TLE1_COL_SATNUM, TLE1_LEN_SATNUM);
    fields_[FLD_INTLDESC].first = line_one_.substr(TLE1_COL_INTLDESC_A, TLE1_LEN_INTLDESC_A +
            TLE1_LEN_INTLDESC_B + TLE1_LEN_INTLDESC_C);
    fields_[FLD_EPOCHYEAR].first = line_one_.substr(TLE1_COL_EPOCH_A, TLE1_LEN_EPOCH_A);

    fields_[FLD_EPOCHDAY].first = line_one_.substr(TLE1_COL_EPOCH_B, TLE1_LEN_EPOCH_B);

    if (line_one_[TLE1_COL_MEANMOTIONDT] == '-') {
        // value is negative
        fields_[FLD_MMOTIONDT].first = "-0";
    } else
        fields_[FLD_MMOTIONDT].first = "0";

    fields_[FLD_MMOTIONDT].first += line_one_.substr(TLE1_COL_MEANMOTIONDT + 1, TLE1_LEN_MEANMOTIONDT);

    // decimal point assumed; exponential notation
    fields_[FLD_MMOTIONDT2].first = ExpToDecimal(line_one_.substr(TLE1_COL_MEANMOTIONDT2, TLE1_LEN_MEANMOTIONDT2));

    // decimal point assumed; exponential notation
    fields_[FLD_BSTAR].first = ExpToDecimal(line_one_.substr(TLE1_COL_BSTAR, TLE1_LEN_BSTAR));

    fields_[FLD_SET].first = line_one_.substr(TLE1_COL_ELNUM, TLE1_LEN_ELNUM);
    TrimLeft(fields_[FLD_SET].first);

    /*
     * line 2
     */

    fields_[FLD_I].first = line_two_.substr(TLE2_COL_INCLINATION, TLE2_LEN_INCLINATION);
    TrimLeft(fields_[FLD_I].first);

    fields_[FLD_RAAN].first = line_two_.substr(TLE2_COL_RAASCENDNODE, TLE2_LEN_RAASCENDNODE);
    TrimLeft(fields_[FLD_RAAN].first);

    // decimal point is assumed
    fields_[FLD_E].first = "0.";
    fields_[FLD_E].first += line_two_.substr(TLE2_COL_ECCENTRICITY, TLE2_LEN_ECCENTRICITY);

    fields_[FLD_ARGPER].first = line_two_.substr(TLE2_COL_ARGPERIGEE, TLE2_LEN_ARGPERIGEE);
    TrimLeft(fields_[FLD_ARGPER].first);

    fields_[FLD_M].first = line_two_.substr(TLE2_COL_MEANANOMALY, TLE2_LEN_MEANANOMALY);
    TrimLeft(fields_[FLD_M].first);

    fields_[FLD_MMOTION].first = line_two_.substr(TLE2_COL_MEANMOTION, TLE2_LEN_MEANMOTION);
    TrimLeft(fields_[FLD_MMOTION].first);

    fields_[FLD_ORBITNUM].first = line_two_.substr(TLE2_COL_REVATEPOCH, TLE2_LEN_REVATEPOCH);
    TrimLeft(fields_[FLD_ORBITNUM].first);

    // update double variables
    for (int fld = FLD_FIRST; fld < FLD_LAST; fld++) {
        fields_[fld].second = atof(fields_[fld].first.c_str());
    }

    int year = static_cast<int> (GetField(Tle::FLD_EPOCHYEAR));
    if (year < 57)
        year += 2000;
    else
        year += 1900;
    double day = GetField(Tle::FLD_EPOCHDAY);

    date_ = Julian(year, day);
}

/*
 * validate a line
 */
bool Tle::IsValidLine(std::string& str, TleLine line) {
    TrimLeft(str);
    TrimRight(str);

    const std::string line1_validation = "1 NNNNNC NNNNAAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN";
    const std::string line2_validation = "2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN";

    if (Tle::LINE_ONE == line) {
        if (!ValidateLine(line1_validation, str))
            return false;
    } else if (Tle::LINE_TWO == line) {
        if (!ValidateLine(line2_validation, str))
            return false;
    }

    // Last char in string must be checksum
    // int chk = CheckSum(str);
    // if (chk != (str[TLE_LEN_LINE_DATA - 1] - '0'))
    //   return false;

    return true;
}

/*
 * validate line given a pattern
 */
bool Tle::ValidateLine(const std::string& pattern, const std::string& line) {

    if (pattern.length() != line.length()) {
        return false;
    }

    std::string::size_type pos = 0;

    while (pos != line.length() - 1) {
        if (isdigit(pattern[pos]) || pattern[pos] == ' ' || pattern[pos] == '.') {
            if (pattern[pos] != line[pos])
                return false;
        } else if (pattern[pos] == 'N') {
            if (!isdigit(line[pos]) && line[pos] != ' ')
                return false;
        }
        pos++;
    }

    return true;
}

/*
 * compute checksum
 */
int Tle::CheckSum(const std::string& str) {
    size_t len = str.size() - 1;
    int xsum = 0;

    for (size_t i = 0; i < len; i++) {
        char ch = str[i];
        if (isdigit(ch))
            xsum += (ch - '0');
        else if (ch == '-')
            xsum++;
    }

    return (xsum % 10);
}

/*
 * trim left of string
 */
void Tle::TrimLeft(std::string& s) {
    if (!s.empty()) {
        std::string::size_type pos = s.find_first_not_of(' ');

        if (pos != std::string::npos)
            s.erase(0, pos);
        else
            s.clear();
    }
}

/*
 * trim right of string
 */
void Tle::TrimRight(std::string& s) {
    if (!s.empty()) {
        std::string::size_type pos = s.find_last_not_of(' ');

        if (pos != std::string::npos)
            s.erase(pos + 1);
        else
            s.clear();
    }
}

