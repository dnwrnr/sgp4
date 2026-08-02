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

#include <cstdio>
#include <sstream>
#include <vector>
#include <locale>

namespace libsgp4
{
namespace
{
    static const unsigned int TLE1_COL_NORADNUM = 2;
    static const unsigned int TLE1_LEN_NORADNUM = 5;
    static const unsigned int TLE1_COL_INTLDESC_A = 9;
    static const unsigned int TLE1_LEN_INTLDESC_A = 2;
//  static const unsigned int TLE1_COL_INTLDESC_B = 11;
    static const unsigned int TLE1_LEN_INTLDESC_B = 3;
//  static const unsigned int TLE1_COL_INTLDESC_C = 14;
    static const unsigned int TLE1_LEN_INTLDESC_C = 3;
    static const unsigned int TLE1_COL_EPOCH_A = 18;
    static const unsigned int TLE1_LEN_EPOCH_A = 2;
    static const unsigned int TLE1_COL_EPOCH_B = 20;
    static const unsigned int TLE1_LEN_EPOCH_B = 12;
    static const unsigned int TLE1_COL_MEANMOTIONDT2 = 33;
    static const unsigned int TLE1_LEN_MEANMOTIONDT2 = 10;
    static const unsigned int TLE1_COL_MEANMOTIONDDT6 = 44;
    static const unsigned int TLE1_LEN_MEANMOTIONDDT6 = 8;
    static const unsigned int TLE1_COL_BSTAR = 53;
    static const unsigned int TLE1_LEN_BSTAR = 8;
//  static const unsigned int TLE1_COL_EPHEMTYPE = 62;
//  static const unsigned int TLE1_LEN_EPHEMTYPE = 1;
//  static const unsigned int TLE1_COL_ELNUM = 64;
//  static const unsigned int TLE1_LEN_ELNUM = 4;

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

/**
 * Initialise the tle object.
 * @exception TleException
 */
void Tle::Initialize()
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

    unsigned int satNumber1;
    unsigned int satNumber2;

    ExtractInteger(mLineOne.substr(TLE1_COL_NORADNUM,
                TLE1_LEN_NORADNUM), satNumber1);
    ExtractInteger(mLineTwo.substr(TLE2_COL_NORADNUM,
                TLE2_LEN_NORADNUM), satNumber2);

    if (satNumber1 != satNumber2)
    {
        throw TleException("Satellite numbers do not match");
    }

    mNoradNumber = satNumber1;

    if (mName.empty())
    {
        mName = mLineOne.substr(TLE1_COL_NORADNUM, TLE1_LEN_NORADNUM);
    }

    mIntDesignator = mLineOne.substr(TLE1_COL_INTLDESC_A,
            TLE1_LEN_INTLDESC_A + TLE1_LEN_INTLDESC_B + TLE1_LEN_INTLDESC_C);

    unsigned int year = 0;
    double day = 0.0;

    ExtractInteger(mLineOne.substr(TLE1_COL_EPOCH_A,
                TLE1_LEN_EPOCH_A), year);
    ExtractDouble(mLineOne.substr(TLE1_COL_EPOCH_B,
                TLE1_LEN_EPOCH_B), 4, day);
    ExtractDouble(mLineOne.substr(TLE1_COL_MEANMOTIONDT2,
                TLE1_LEN_MEANMOTIONDT2), 2, mMeanMotionDt2);
    ExtractExponential(mLineOne.substr(TLE1_COL_MEANMOTIONDDT6,
                TLE1_LEN_MEANMOTIONDDT6), mMeanMotionDdt6);
    ExtractExponential(mLineOne.substr(TLE1_COL_BSTAR,
                TLE1_LEN_BSTAR), mBstar);

    /*
     * line 2
     */
    ExtractDouble(mLineTwo.substr(TLE2_COL_INCLINATION,
                TLE2_LEN_INCLINATION), 4, mInclination);
    ExtractDouble(mLineTwo.substr(TLE2_COL_RAASCENDNODE,
                TLE2_LEN_RAASCENDNODE), 4, mRightAscendingNode);
    ExtractDouble(mLineTwo.substr(TLE2_COL_ECCENTRICITY,
                TLE2_LEN_ECCENTRICITY), -1, mEccentricity);
    ExtractDouble(mLineTwo.substr(TLE2_COL_ARGPERIGEE,
                TLE2_LEN_ARGPERIGEE), 4, mArgumentPerigee);
    ExtractDouble(mLineTwo.substr(TLE2_COL_MEANANOMALY,
                TLE2_LEN_MEANANOMALY), 4, mMeanAnomaly);
    ExtractDouble(mLineTwo.substr(TLE2_COL_MEANMOTION,
                TLE2_LEN_MEANMOTION), 3, mMeanMotion);
    ExtractInteger(mLineTwo.substr(TLE2_COL_REVATEPOCH,
                TLE2_LEN_REVATEPOCH), mOrbitNumber);
    
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
    return str.length() == LineLength() ? true : false;
}

/**
 * Convert a string containing an integer
 * @param[in] str The string to convert
 * @param[out] val The result
 * @exception TleException on conversion error
 */
void Tle::ExtractInteger(const std::string& str, unsigned int& val)
{
    bool foundDigit = false;
    unsigned int temp = 0;

    for (auto& i : str)
    {
        if (isdigit(static_cast<unsigned char>(i)))
        {
            foundDigit = true;
            temp = (temp * 10) + static_cast<unsigned int>(i - '0');
        }
        else if (foundDigit)
        {
            throw TleException("Unexpected non digit");
        }
        else if (i != ' ')
        {
            throw TleException("Invalid character");
        }
    }

    if (!foundDigit)
    {
        val = 0;
    }
    else
    {
        val = temp;
    }
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
                if(*i == '-' || *i == '+')
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
                 */
                temp += '0';
                temp += '.';
            }
            
            /*
             * should be a digit
             */
            if (isdigit(*i))
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
            if (isdigit(*i))
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
}

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
    double meanMotion = std::stod(fields[3]);
    double eccentricity = std::stod(fields[4]);
    double inclination = std::stod(fields[5]);
    double raan = std::stod(fields[6]);
    double argPerigee = std::stod(fields[7]);
    double meanAnomaly = std::stod(fields[8]);
    // fields[9] = ephemeris type (unused)
    // fields[10] = classification type (unused)
    unsigned int noradNumber = static_cast<unsigned int>(std::stoul(fields[11]));
    // fields[12] = element set number (unused)
    unsigned int orbitNumber = static_cast<unsigned int>(std::stoul(fields[13]));
    double bstar = std::stod(fields[14]);
    double meanMotionDt2 = std::stod(fields[15]);
    double meanMotionDdt6 = std::stod(fields[16]);

    int year = 0, month = 0, day = 0, hour = 0, minute = 0, second = 0;
    if (std::sscanf(epochStr.c_str(), "%d-%d-%dT%d:%d:%d",
                     &year, &month, &day, &hour, &minute, &second) != 6)
    {
        throw TleException("Invalid epoch format");
    }
    std::string::size_type dotPos = epochStr.rfind('.');
    int microsecond = dotPos == std::string::npos ? 0 :
        ParseIsoMicrosecond(epochStr.substr(dotPos + 1));

    DateTime epoch(year, month, day, hour, minute, second, microsecond);

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
