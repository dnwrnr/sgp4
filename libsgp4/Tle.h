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


#pragma once

#include "DateTime.h"
#include "TleException.h"
#include "Util.h"

namespace libsgp4
{

/**
 * @brief Processes a two-line element set used to convey OrbitalElements.
 *
 * Used to extract the various raw fields from a two-line element set.
 */
class Tle
{
public:
    /**
     * @details Initialise given the two lines of a tle
     * @param[in] lineOne Tle line one
     * @param[in] lineTwo Tle line two
     */
    Tle(std::string lineOne, std::string lineTwo)
        : mLineOne(std::move(lineOne))
        , mLineTwo(std::move(lineTwo))
    {
        Initialize();
    }

    /**
     * @details Initialise given the satellite name and the two lines of a tle
     * @param[in] name Satellite name
     * @param[in] lineOne Tle line one
     * @param[in] lineTwo Tle line two
     */
    Tle(std::string name, std::string lineOne, std::string lineTwo)
        : mName(std::move(name))
        , mLineOne(std::move(lineOne))
        , mLineTwo(std::move(lineTwo))
    {
        Initialize();
    }

    /**
     * @details Construct a Tle from a CelesTrak CSV data line.
     *
     * Expected CSV columns (header row is not parsed here):
     *   OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,
     *   INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,
     *   EPHEMERIS_TYPE,CLASSIFICATION_TYPE,NORAD_CAT_ID,ELEMENT_SET_NO,
     *   REV_AT_EPOCH,BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT
     *
     * @param csvLine A single CSV data line (no header)
     * @returns A Tle object with all orbital elements populated.
     *          Line1() and Line2() will return empty strings.
     * @throws TleException if the line is malformed or has wrong field count
     */
    static Tle FromCsv(const std::string& csvLine);

    /**
     * Copy constructor
     * @param[in] tle Tle object to copy from
     */
    Tle(const Tle& tle)
    {
        mName = tle.mName;
        mLineOne = tle.mLineOne;
        mLineTwo = tle.mLineTwo;

        mNoradNumber = tle.mNoradNumber;
        mIntDesignator = tle.mIntDesignator;
        mEpoch = tle.mEpoch;
        mMeanMotionDt2 = tle.mMeanMotionDt2;
        mMeanMotionDdt6 = tle.mMeanMotionDdt6;
        mBstar = tle.mBstar;
        mInclination = tle.mInclination;
        mRightAscendingNode = tle.mRightAscendingNode;
        mEccentricity = tle.mEccentricity;
        mArgumentPerigee = tle.mArgumentPerigee;
        mMeanAnomaly = tle.mMeanAnomaly;
        mMeanMotion = tle.mMeanMotion;
        mOrbitNumber = tle.mOrbitNumber;
    }

    /**
     * Get the satellite name
     * @returns the satellite name
     */
    std::string Name() const
    {
        return mName;
    }

    /**
     * Get the first line of the tle
     * @returns the first line of the tle
     */
    std::string Line1() const
    {
        return mLineOne;
    }

    /**
     * Get the second line of the tle
     * @returns the second line of the tle
     */
    std::string Line2() const
    {
        return mLineTwo;
    }

    /**
     * Get the norad number
     * @returns the norad number
     */
    unsigned int NoradNumber() const
    {
        return mNoradNumber;
    }

    /**
     * Get the international designator
     * @returns the international designator
     */
    std::string IntDesignator() const
    {
        return mIntDesignator;
    }

    /**
     * Get the tle epoch
     * @returns the tle epoch
     */
    DateTime Epoch() const
    {
        return mEpoch;
    }

    /**
     * Get the first time derivative of the mean motion divided by two
     * @returns the first time derivative of the mean motion divided by two
     */
    double MeanMotionDt2() const
    {
        return mMeanMotionDt2;
    }

    /**
     * Get the second time derivative of mean motion divided by six
     * @returns the second time derivative of mean motion divided by six
     */
    double MeanMotionDdt6() const
    {
        return mMeanMotionDdt6;
    }

    /**
     * Get the BSTAR drag term
     * @returns the BSTAR drag term
     */
    double BStar() const
    {
        return mBstar;
    }

    /**
     * Get the inclination
     * @param inDegrees Whether to return the value in degrees or radians
     * @returns the inclination
     */
    double Inclination(bool inDegrees) const
    {
        if (inDegrees)
        {
            return mInclination;
        }
        else
        {
            return Util::DegreesToRadians(mInclination);
        }
    }

    /**
     * Get the right ascension of the ascending node
     * @param inDegrees Whether to return the value in degrees or radians
     * @returns the right ascension of the ascending node
     */
    double RightAscendingNode(const bool inDegrees) const
    {
        if (inDegrees)
        {
            return mRightAscendingNode;
        }
        else
        {
            return Util::DegreesToRadians(mRightAscendingNode);
        }
    }

    /**
     * Get the eccentricity
     * @returns the eccentricity
     */
    double Eccentricity() const
    {
        return mEccentricity;
    }

    /**
     * Get the argument of perigee
     * @param inDegrees Whether to return the value in degrees or radians
     * @returns the argument of perigee
     */
    double ArgumentPerigee(const bool inDegrees) const
    {
        if (inDegrees)
        {
            return mArgumentPerigee;
        }
        else
        {
            return Util::DegreesToRadians(mArgumentPerigee);
        }
    }

    /**
     * Get the mean anomaly
     * @param inDegrees Whether to return the value in degrees or radians
     * @returns the mean anomaly
     */
    double MeanAnomaly(const bool inDegrees) const
    {
        if (inDegrees)
        {
            return mMeanAnomaly;
        }
        else
        {
            return Util::DegreesToRadians(mMeanAnomaly);
        }
    }

    /**
     * Get the mean motion
     * @returns the mean motion (revolutions per day)
     */
    double MeanMotion() const
    {
        return mMeanMotion;
    }

    /**
     * Get the orbit number
     * @returns the orbit number
     */
    unsigned int OrbitNumber() const
    {
        return mOrbitNumber;
    }

    /**
     * Get the expected tle line length
     * @returns the tle line length
     */
    static unsigned int LineLength()
    {
        return TLE_LEN_LINE_DATA;
    }
    
    /**
     * Dump this object to a string
     * @returns string
     */
    std::string ToString() const
    {
        std::stringstream ss;
        ss << std::right << std::fixed;
        ss << "Norad Number:         " << NoradNumber() << std::endl;
        ss << "Int. Designator:      " << IntDesignator() << std::endl;
        ss << "Epoch:                " << Epoch() << std::endl;
        ss << "Orbit Number:         " << OrbitNumber() << std::endl;
        ss << std::setprecision(8);
        ss << "Mean Motion Dt2:      ";
        ss << std::setw(12) << MeanMotionDt2() << std::endl;
        ss << "Mean Motion Ddt6:     ";
        ss << std::setw(12) << MeanMotionDdt6() << std::endl;
        ss << "Eccentricity:         ";
        ss << std::setw(12) << Eccentricity() << std::endl;
        ss << "BStar:                ";
        ss << std::setw(12) << BStar() << std::endl;
        ss << "Inclination:          ";
        ss << std::setw(12) << Inclination(true) << std::endl;
        ss << "Right Ascending Node: ";
        ss << std::setw(12) << RightAscendingNode(true) << std::endl;
        ss << "Argument Perigee:     ";
        ss << std::setw(12) << ArgumentPerigee(true) << std::endl;
        ss << "Mean Anomaly:         ";
        ss << std::setw(12) << MeanAnomaly(true) << std::endl;
        ss << "Mean Motion:          ";
        ss << std::setw(12) << MeanMotion() << std::endl;
        return ss.str();
    }

private:
    Tle(const std::string& name,
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
        unsigned int orbitNumber);

    void Initialize();
    static bool IsValidLineLength(const std::string& str);
    void ExtractInteger(const std::string& str, unsigned int& val);
    void ExtractDouble(const std::string& str, int pointPos, double& val);
    void ExtractExponential(const std::string& str, double& val);

private:
    std::string mName;
    std::string mLineOne;
    std::string mLineTwo;

    std::string mIntDesignator;
    DateTime mEpoch;
    double mMeanMotionDt2{};
    double mMeanMotionDdt6{};
    double mBstar{};
    double mInclination{};
    double mRightAscendingNode{};
    double mEccentricity{};
    double mArgumentPerigee{};
    double mMeanAnomaly{};
    double mMeanMotion{};
    unsigned int mNoradNumber{};
    unsigned int mOrbitNumber{};

    static const unsigned int TLE_LEN_LINE_DATA = 69;
    static const unsigned int TLE_LEN_LINE_NAME = 22;
};


inline std::ostream& operator<<(std::ostream& strm, const Tle& t)
{
    return strm << t.ToString();
}

} // namespace libsgp4
