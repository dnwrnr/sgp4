/*
 * Copyright 2013 Daniel Warner <contact@danrw.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "Tle.h"

#include <charconv>
#include <cmath>
#include <locale>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace libsgp4
{
namespace
{
    const unsigned int TLE1_COL_NORADNUM = 2;
    const unsigned int TLE1_LEN_NORADNUM = 5;
    const unsigned int TLE1_COL_INTLDESC_A = 9;
    const unsigned int TLE1_LEN_INTLDESC_A = 2;
    //  static const unsigned int TLE1_COL_INTLDESC_B = 11;
    const unsigned int TLE1_LEN_INTLDESC_B = 3;
    //  static const unsigned int TLE1_COL_INTLDESC_C = 14;
    const unsigned int TLE1_LEN_INTLDESC_C = 3;
    const unsigned int TLE1_COL_EPOCH_A = 18;
    const unsigned int TLE1_LEN_EPOCH_A = 2;
    const unsigned int TLE1_COL_EPOCH_B = 20;
    const unsigned int TLE1_LEN_EPOCH_B = 12;
    const unsigned int TLE1_COL_MEANMOTIONDT2 = 33;
    const unsigned int TLE1_LEN_MEANMOTIONDT2 = 10;
    const unsigned int TLE1_COL_MEANMOTIONDDT6 = 44;
    const unsigned int TLE1_LEN_MEANMOTIONDDT6 = 8;
    const unsigned int TLE1_COL_BSTAR = 53;
    const unsigned int TLE1_LEN_BSTAR = 8;
    //  static const unsigned int TLE1_COL_EPHEMTYPE = 62;
    //  static const unsigned int TLE1_LEN_EPHEMTYPE = 1;
    //  static const unsigned int TLE1_COL_ELNUM = 64;
    //  static const unsigned int TLE1_LEN_ELNUM = 4;

    const unsigned int TLE2_COL_NORADNUM = 2;
    const unsigned int TLE2_LEN_NORADNUM = 5;
    const unsigned int TLE2_COL_INCLINATION = 8;
    const unsigned int TLE2_LEN_INCLINATION = 8;
    const unsigned int TLE2_COL_RAASCENDNODE = 17;
    const unsigned int TLE2_LEN_RAASCENDNODE = 8;
    const unsigned int TLE2_COL_ECCENTRICITY = 26;
    const unsigned int TLE2_LEN_ECCENTRICITY = 7;
    const unsigned int TLE2_COL_ARGPERIGEE = 34;
    const unsigned int TLE2_LEN_ARGPERIGEE = 8;
    const unsigned int TLE2_COL_MEANANOMALY = 43;
    const unsigned int TLE2_LEN_MEANANOMALY = 8;
    const unsigned int TLE2_COL_MEANMOTION = 52;
    const unsigned int TLE2_LEN_MEANMOTION = 11;
    const unsigned int TLE2_COL_REVATEPOCH = 63;
    const unsigned int TLE2_LEN_REVATEPOCH = 5;

    // Alpha-5 prefixes ordered by ascending value. The letters I and O are omitted to avoid
    // confusion with the digits 1 and 0, so the first entry maps to a leading value of 10.
    const char* const ALPHA5_PREFIXES = "ABCDEFGHJKLMNPQRSTUVWXYZ";
    const unsigned int ALPHA5_FIRST_LEADING_VALUE = 10;
    const unsigned int ALPHA5_TAIL_SCALE = 10000;

    /**
     * Convert a field of decimal digits to an unsigned integer.
     *
     * Unlike std::stoul this neither skips leading whitespace nor accepts a sign, and it
     * reports an out of range field instead of wrapping around, so the caller sees exactly
     * the digits that were in the field.
     *
     * @param[in] text The digits to convert
     * @param[in] description The field name to report on failure
     * @returns The converted value
     * @exception TleException if the field is empty, holds a non digit or does not fit
     */
    unsigned int ParseDigits(std::string_view text, const char* description)
    {
        unsigned int value = 0;
        const std::from_chars_result result = std::from_chars(text.data(), text.data() + text.size(), value);
        if (result.ec != std::errc() || result.ptr != text.data() + text.size())
        {
            throw TleException(description);
        }
        return value;
    }

    /**
     * Decode an Alpha-5 object number from the satellite number field of a tle.
     *
     * Alpha-5 replaces the leading digit of object numbers from 100000 upwards with a letter,
     * so A0000 is 100000 and Z9999 is 339999. Object numbers below 100000 are unaffected.
     *
     * @param[in] field The satellite number field
     * @param[out] val The decoded object number
     * @returns Whether the field held an Alpha-5 object number
     * @exception TleException on an unsupported prefix or a non digit tail
     */
    bool DecodeAlpha5NoradNumber(const std::string& field, unsigned int& val)
    {
        if (field.empty() || field[0] < 'A' || field[0] > 'Z')
        {
            return false;
        }

        const std::string_view prefixes(ALPHA5_PREFIXES);
        const std::string_view::size_type index = prefixes.find(field[0]);

        if (index == std::string_view::npos)
        {
            throw TleException("Unsupported Alpha-5 satellite number prefix");
        }

        if (field.length() != TLE1_LEN_NORADNUM)
        {
            throw TleException("Invalid length for Alpha-5 satellite number");
        }

        const unsigned int tail = ParseDigits(std::string_view(field).substr(1), "Invalid Alpha-5 satellite number");

        const unsigned int leading = ALPHA5_FIRST_LEADING_VALUE + static_cast<unsigned int>(index);
        val = (leading * ALPHA5_TAIL_SCALE) + tail;

        return true;
    }
} // namespace

/**
 * Initialise the tle object.
 * @exception TleException
 */
void Tle::Initialise()
{
    if (!IsValidLineLength(mLineOne))
    {
        throw TleException("Invalid length for line one");
    }

    if (!IsValidLineLength(mLineTwo))
    {
        throw TleException("Invalid length for line two");
    }

    if (mLineOne[0] != '1')
    {
        throw TleException("Invalid line beginning for line one");
    }

    if (mLineTwo[0] != '2')
    {
        throw TleException("Invalid line beginning for line two");
    }

    unsigned int satNumber1 = 0;
    unsigned int satNumber2 = 0;

    const std::string satNumberField1 = mLineOne.substr(TLE1_COL_NORADNUM, TLE1_LEN_NORADNUM);
    const std::string satNumberField2 = mLineTwo.substr(TLE2_COL_NORADNUM, TLE2_LEN_NORADNUM);

    if (!DecodeAlpha5NoradNumber(satNumberField1, satNumber1))
    {
        ExtractInteger(satNumberField1, satNumber1);
    }

    if (!DecodeAlpha5NoradNumber(satNumberField2, satNumber2))
    {
        ExtractInteger(satNumberField2, satNumber2);
    }

    if (satNumber1 != satNumber2)
    {
        throw TleException("Satellite numbers do not match");
    }

    mNoradNumber = satNumber1;

    if (mName.empty())
    {
        mName = mLineOne.substr(TLE1_COL_NORADNUM, TLE1_LEN_NORADNUM);
    }

    mIntDesignator =
        mLineOne.substr(TLE1_COL_INTLDESC_A, TLE1_LEN_INTLDESC_A + TLE1_LEN_INTLDESC_B + TLE1_LEN_INTLDESC_C);

    unsigned int year = 0;
    double day = 0.0;

    ExtractInteger(mLineOne.substr(TLE1_COL_EPOCH_A, TLE1_LEN_EPOCH_A), year);
    ExtractDouble(mLineOne.substr(TLE1_COL_EPOCH_B, TLE1_LEN_EPOCH_B), 4, day);
    ExtractDouble(mLineOne.substr(TLE1_COL_MEANMOTIONDT2, TLE1_LEN_MEANMOTIONDT2), 2, mMeanMotionDt2);
    ExtractExponential(mLineOne.substr(TLE1_COL_MEANMOTIONDDT6, TLE1_LEN_MEANMOTIONDDT6), mMeanMotionDdt6);
    ExtractExponential(mLineOne.substr(TLE1_COL_BSTAR, TLE1_LEN_BSTAR), mBstar);

    /*
     * line 2
     */
    ExtractDouble(mLineTwo.substr(TLE2_COL_INCLINATION, TLE2_LEN_INCLINATION), 4, mInclination);
    ExtractDouble(mLineTwo.substr(TLE2_COL_RAASCENDNODE, TLE2_LEN_RAASCENDNODE), 4, mRightAscendingNode);
    ExtractDouble(mLineTwo.substr(TLE2_COL_ECCENTRICITY, TLE2_LEN_ECCENTRICITY), -1, mEccentricity);
    ExtractDouble(mLineTwo.substr(TLE2_COL_ARGPERIGEE, TLE2_LEN_ARGPERIGEE), 4, mArgumentPerigee);
    ExtractDouble(mLineTwo.substr(TLE2_COL_MEANANOMALY, TLE2_LEN_MEANANOMALY), 4, mMeanAnomaly);
    ExtractDouble(mLineTwo.substr(TLE2_COL_MEANMOTION, TLE2_LEN_MEANMOTION), 3, mMeanMotion);
    ExtractInteger(mLineTwo.substr(TLE2_COL_REVATEPOCH, TLE2_LEN_REVATEPOCH), mOrbitNumber);

    if (year < 57)
    {
        year += 2000;
    }
    else
    {
        year += 1900;
    }

    mEpoch = DateTime(year, day);
}

/**
 * Check
 * @param str The string to check
 * @returns Whether true of the string has a valid length
 */
bool Tle::IsValidLineLength(const std::string& str)
{
    return str.size() == LineLength();
}

/**
 * Convert a string containing an integer
 *
 * The field is either blank, in which case it is zero, or space padded digits. A space after
 * the first digit is a misplaced field rather than padding, so it is rejected.
 *
 * @param[in] str The string to convert
 * @param[out] val The result
 * @exception TleException on conversion error
 */
void Tle::ExtractInteger(const std::string& str, unsigned int& val)
{
    const std::string_view field(str);
    const std::string_view::size_type firstDigit = field.find_first_not_of(' ');

    if (firstDigit == std::string_view::npos)
    {
        val = 0;
        return;
    }

    val = ParseDigits(field.substr(firstDigit), "Unexpected non digit");
}

/**
 * Convert a string containing an double
 * @param[in] str The string to convert
 * @param[in] pointPos The position of the decimal point. (-1 if none)
 * @param[out] val The result
 * @exception TleException on conversion error
 */
void Tle::ExtractDouble(const std::string& str, int pointPos, double& val)
{
    std::string temp;
    bool foundDigit = false;

    for (std::string::const_iterator i = str.begin(); i != str.end(); ++i)
    {
        /*
         * integer part
         */
        if (pointPos >= 0 && i < str.begin() + pointPos - 1)
        {
            bool done = false;

            if (i == str.begin())
            {
                if (*i == '-' || *i == '+')
                {
                    /*
                     * first character could be signed
                     */
                    temp += *i;
                    done = true;
                }
            }

            if (!done)
            {
                if (isdigit(static_cast<unsigned char>(*i)))
                {
                    foundDigit = true;
                    temp += *i;
                }
                else if (foundDigit)
                {
                    throw TleException("Unexpected non digit");
                }
                else if (*i != ' ')
                {
                    throw TleException("Invalid character");
                }
            }
        }
        /*
         * decimal point
         */
        else if (pointPos >= 0 && i == str.begin() + pointPos - 1)
        {
            if (temp.length() == 0)
            {
                /*
                 * integer part is blank, so add a '0'
                 */
                temp += '0';
            }

            if (*i == '.')
            {
                /*
                 * decimal point found
                 */
                temp += *i;
            }
            else
            {
                throw TleException("Failed to find decimal point");
            }
        }
        /*
         * fraction part
         */
        else
        {
            if (i == str.begin() && pointPos == -1)
            {
                /*
                 * no decimal point expected, add 0. beginning
                 *
                 * this is the eccentricity field of line 2, seven digits with a leading
                 * decimal point assumed. an eccentricity always fits in seven digits, so the
                 * field carries no padding and a space here is a malformed line rather than
                 * a blank leading digit.
                 */
                temp += '0';
                temp += '.';
            }

            /*
             * should be a digit
             */
            if (isdigit(static_cast<unsigned char>(*i)))
            {
                temp += *i;
            }
            else
            {
                throw TleException("Invalid digit");
            }
        }
    }

    if (!Util::FromString<double>(temp, val))
    {
        throw TleException("Failed to convert value to double");
    }
}

/**
 * Convert a string containing an exponential
 * @param[in] str The string to convert
 * @param[out] val The result
 * @exception TleException on conversion error
 */
void Tle::ExtractExponential(const std::string& str, double& val)
{
    std::string temp;

    for (std::string::const_iterator i = str.begin(); i != str.end(); ++i)
    {
        if (i == str.begin())
        {
            if (*i == '-' || *i == '+' || *i == ' ')
            {
                if (*i == '-')
                {
                    temp += *i;
                }
                temp += '0';
                temp += '.';
            }
            else
            {
                throw TleException("Invalid sign");
            }
        }
        else if (i == str.end() - 2)
        {
            if (*i == '-' || *i == '+')
            {
                temp += 'e';
                temp += *i;
            }
            else
            {
                throw TleException("Invalid exponential sign");
            }
        }
        else
        {
            if (isdigit(static_cast<unsigned char>(*i)))
            {
                temp += *i;
            }
            else
            {
                throw TleException("Invalid digit");
            }
        }
    }

    if (!Util::FromString<double>(temp, val))
    {
        throw TleException("Failed to convert value to double");
    }
}

/**
 * Construct a Tle directly from parsed fields.
 */
Tle::Tle(const std::string& name,
         unsigned int noradNumber,
         const std::string& intDesignator,
         const DateTime& epoch,
         double meanMotionDt2,
         double meanMotionDdt6,
         double bstar,
         double inclination,
         double rightAscendingNode,
         double eccentricity,
         double argumentPerigee,
         double meanAnomaly,
         double meanMotion,
         unsigned int orbitNumber)
    : mName(name)
    , mIntDesignator(intDesignator)
    , mEpoch(epoch)
    , mMeanMotionDt2(meanMotionDt2)
    , mMeanMotionDdt6(meanMotionDdt6)
    , mBstar(bstar)
    , mInclination(inclination)
    , mRightAscendingNode(rightAscendingNode)
    , mEccentricity(eccentricity)
    , mArgumentPerigee(argumentPerigee)
    , mMeanAnomaly(meanAnomaly)
    , mMeanMotion(meanMotion)
    , mNoradNumber(noradNumber)
    , mOrbitNumber(orbitNumber)
{
}

namespace
{
    std::vector<std::string> SplitCsv(const std::string& line)
    {
        std::vector<std::string> fields;
        std::istringstream stream(line);
        std::string field;
        while (std::getline(stream, field, ','))
        {
            if (!field.empty() && field.back() == '\r')
            {
                field.pop_back();
            }
            fields.push_back(std::move(field));
        }
        return fields;
    }

    int ParseIsoMicrosecond(const std::string& s)
    {
        if (s.empty())
        {
            return 0;
        }
        std::string digits;
        for (char c : s)
        {
            if (std::isdigit(static_cast<unsigned char>(c)))
            {
                digits += c;
            }
        }
        if (digits.empty())
        {
            return 0;
        }
        while (digits.length() < 6)
        {
            digits += '0';
        }
        return std::stoi(digits.substr(0, 6));
    }

    bool ParseUnsigned(const std::string& s, std::string::size_type& pos, unsigned int& value)
    {
        const char* const begin = s.data() + pos;
        const char* const end = s.data() + s.size();
        const std::from_chars_result result = std::from_chars(begin, end, value);
        if (result.ec != std::errc() || result.ptr == begin)
        {
            return false;
        }
        pos = static_cast<std::string::size_type>(result.ptr - s.data());
        return true;
    }

    bool ParseIsoDateTime(const std::string& s,
                          unsigned int& year,
                          unsigned int& month,
                          unsigned int& day,
                          unsigned int& hour,
                          unsigned int& minute,
                          unsigned int& second)
    {
        static const char DELIMITERS[] = {'-', '-', 'T', ':', ':'};
        unsigned int* const fields[] = {&year, &month, &day, &hour, &minute, &second};

        std::string::size_type pos = 0;
        for (size_t i = 0; i < 6; ++i)
        {
            if (!ParseUnsigned(s, pos, *fields[i]))
            {
                return false;
            }
            if (i + 1 < 6)
            {
                if (pos >= s.size() || s[pos] != DELIMITERS[i])
                {
                    return false;
                }
                ++pos;
            }
        }
        return pos == s.size();
    }

    double ParseCsvDouble(const std::string& field, const char* description)
    {
        std::string::size_type pos = 0;
        double value = 0.0;
        try
        {
            value = std::stod(field, &pos);
        }
        catch (const std::invalid_argument&)
        {
            throw TleException((std::string("Invalid CSV value for ") + description).c_str());
        }
        catch (const std::out_of_range&)
        {
            throw TleException((std::string("CSV value out of range for ") + description).c_str());
        }
        if (pos != field.size())
        {
            throw TleException((std::string("Invalid trailing characters in CSV ") + description).c_str());
        }
        return value;
    }

    unsigned int ParseCsvUnsigned(const std::string& field, const char* description)
    {
        std::string::size_type pos = 0;
        unsigned long value = 0;
        try
        {
            value = std::stoul(field, &pos);
        }
        catch (const std::invalid_argument&)
        {
            throw TleException((std::string("Invalid CSV value for ") + description).c_str());
        }
        catch (const std::out_of_range&)
        {
            throw TleException((std::string("CSV value out of range for ") + description).c_str());
        }
        if (pos != field.size())
        {
            throw TleException((std::string("Invalid trailing characters in CSV ") + description).c_str());
        }
        return static_cast<unsigned int>(value);
    }
} // namespace

Tle Tle::FromCsv(const std::string& csvLine)
{
    const unsigned int EXPECTED_FIELDS = 17;
    std::vector<std::string> fields = SplitCsv(csvLine);

    if (fields.size() != EXPECTED_FIELDS)
    {
        throw TleException("Invalid CSV field count");
    }

    const std::string& name = fields[0];
    const std::string& intDesignator = fields[1];
    const std::string& epochStr = fields[2];
    const double meanMotion = ParseCsvDouble(fields[3], "mean motion");
    const double eccentricity = ParseCsvDouble(fields[4], "eccentricity");
    const double inclination = ParseCsvDouble(fields[5], "inclination");
    const double raan = ParseCsvDouble(fields[6], "right ascension");
    const double argPerigee = ParseCsvDouble(fields[7], "argument of perigee");
    const double meanAnomaly = ParseCsvDouble(fields[8], "mean anomaly");
    // fields[9] = ephemeris type (unused)
    // fields[10] = classification type (unused)
    const unsigned int noradNumber = ParseCsvUnsigned(fields[11], "norad number");
    // fields[12] = element set number (unused)
    const unsigned int orbitNumber = ParseCsvUnsigned(fields[13], "orbit number");
    const double bstar = ParseCsvDouble(fields[14], "bstar");
    const double meanMotionDt2 = ParseCsvDouble(fields[15], "mean motion dt2");
    const double meanMotionDdt6 = ParseCsvDouble(fields[16], "mean motion ddt6");

    if (!std::isfinite(meanMotion) || meanMotion <= 0.0)
    {
        throw TleException("Invalid CSV mean motion");
    }
    if (!std::isfinite(eccentricity) || eccentricity < 0.0 || eccentricity >= 1.0)
    {
        throw TleException("Invalid CSV eccentricity");
    }
    if (!std::isfinite(inclination) || inclination < 0.0 || inclination > 180.0)
    {
        throw TleException("Invalid CSV inclination");
    }
    if (!std::isfinite(raan) || !std::isfinite(argPerigee) || !std::isfinite(meanAnomaly))
    {
        throw TleException("Invalid CSV angle");
    }
    if (!std::isfinite(bstar) || !std::isfinite(meanMotionDt2) || !std::isfinite(meanMotionDdt6))
    {
        throw TleException("Invalid CSV drag coefficient");
    }

    std::string::size_type dotPos = epochStr.rfind('.');
    const std::string dateTimeStr = dotPos == std::string::npos ? epochStr : epochStr.substr(0, dotPos);

    unsigned int year = 0, month = 0, day = 0, hour = 0, minute = 0, second = 0;
    if (!ParseIsoDateTime(dateTimeStr, year, month, day, hour, minute, second))
    {
        throw TleException("Invalid epoch format");
    }
    int microsecond = dotPos == std::string::npos ? 0 : ParseIsoMicrosecond(epochStr.substr(dotPos + 1));

    DateTime epoch(static_cast<int>(year),
                   static_cast<int>(month),
                   static_cast<int>(day),
                   static_cast<int>(hour),
                   static_cast<int>(minute),
                   static_cast<int>(second),
                   microsecond);

    return Tle(name,
               noradNumber,
               intDesignator,
               epoch,
               meanMotionDt2,
               meanMotionDdt6,
               bstar,
               inclination,
               raan,
               eccentricity,
               argPerigee,
               meanAnomaly,
               meanMotion,
               orbitNumber);
}

} // namespace libsgp4
