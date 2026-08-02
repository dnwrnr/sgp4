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

#include "Util.h"
#include "DateTime.h"
#include "TleException.h"

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
     * @param[in] line_one Tle line one
     * @param[in] line_two Tle line two
     */
    Tle(std::string line_one, std::string line_two)
        : m_line_one(std::move(line_one))
        , m_line_two(std::move(line_two))
    {
        Initialize();
    }

    /**
     * @details Initialise given the satellite name and the two lines of a tle
     * @param[in] name Satellite name
     * @param[in] line_one Tle line one
     * @param[in] line_two Tle line two
     */
    Tle(std::string name, std::string line_one, std::string line_two)
        : m_name(std::move(name))
        , m_line_one(std::move(line_one))
        , m_line_two(std::move(line_two))
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
     * @param csv_line A single CSV data line (no header)
     * @returns A Tle object with all orbital elements populated.
     *          Line1() and Line2() will return empty strings.
     * @throws TleException if the line is malformed or has wrong field count
     */
    static Tle FromCsv(const std::string& csv_line);

    /**
     * Copy constructor
     * @param[in] tle Tle object to copy from
     */
    Tle(const Tle& tle)
    {
        m_name = tle.m_name;
        m_line_one = tle.m_line_one;
        m_line_two = tle.m_line_two;

        m_norad_number = tle.m_norad_number;
        m_int_designator = tle.m_int_designator;
        m_epoch = tle.m_epoch;
        m_mean_motion_dt2 = tle.m_mean_motion_dt2;
        m_mean_motion_ddt6 = tle.m_mean_motion_ddt6;
        m_bstar = tle.m_bstar;
        m_inclination = tle.m_inclination;
        m_right_ascending_node = tle.m_right_ascending_node;
        m_eccentricity = tle.m_eccentricity;
        m_argument_perigee = tle.m_argument_perigee;
        m_mean_anomaly = tle.m_mean_anomaly;
        m_mean_motion = tle.m_mean_motion;
        m_orbit_number = tle.m_orbit_number;
    }

    /**
     * Get the satellite name
     * @returns the satellite name
     */
    std::string Name() const
    {
        return m_name;
    }

    /**
     * Get the first line of the tle
     * @returns the first line of the tle
     */
    std::string Line1() const
    {
        return m_line_one;
    }

    /**
     * Get the second line of the tle
     * @returns the second line of the tle
     */
    std::string Line2() const
    {
        return m_line_two;
    }

    /**
     * Get the norad number
     * @returns the norad number
     */
    unsigned int NoradNumber() const
    {
        return m_norad_number;
    }

    /**
     * Get the international designator
     * @returns the international designator
     */
    std::string IntDesignator() const
    {
        return m_int_designator;
    }

    /**
     * Get the tle epoch
     * @returns the tle epoch
     */
    DateTime Epoch() const
    {
        return m_epoch;
    }

    /**
     * Get the first time derivative of the mean motion divided by two
     * @returns the first time derivative of the mean motion divided by two
     */
    double MeanMotionDt2() const
    {
        return m_mean_motion_dt2;
    }

    /**
     * Get the second time derivative of mean motion divided by six
     * @returns the second time derivative of mean motion divided by six
     */
    double MeanMotionDdt6() const
    {
        return m_mean_motion_ddt6;
    }

    /**
     * Get the BSTAR drag term
     * @returns the BSTAR drag term
     */
    double BStar() const
    {
        return m_bstar;
    }

    /**
     * Get the inclination
     * @param in_degrees Whether to return the value in degrees or radians
     * @returns the inclination
     */
    double Inclination(bool in_degrees) const
    {
        if (in_degrees)
        {
            return m_inclination;
        }
        else
        {
            return Util::DegreesToRadians(m_inclination);
        }
    }

    /**
     * Get the right ascension of the ascending node
     * @param in_degrees Whether to return the value in degrees or radians
     * @returns the right ascension of the ascending node
     */
    double RightAscendingNode(const bool in_degrees) const
    {
        if (in_degrees)
        {
            return m_right_ascending_node;
        }
        else
        {
            return Util::DegreesToRadians(m_right_ascending_node);
        }
    }

    /**
     * Get the eccentricity
     * @returns the eccentricity
     */
    double Eccentricity() const
    {
        return m_eccentricity;
    }

    /**
     * Get the argument of perigee
     * @param in_degrees Whether to return the value in degrees or radians
     * @returns the argument of perigee
     */
    double ArgumentPerigee(const bool in_degrees) const
    {
        if (in_degrees)
        {
            return m_argument_perigee;
        }
        else
        {
            return Util::DegreesToRadians(m_argument_perigee);
        }
    }

    /**
     * Get the mean anomaly
     * @param in_degrees Whether to return the value in degrees or radians
     * @returns the mean anomaly
     */
    double MeanAnomaly(const bool in_degrees) const
    {
        if (in_degrees)
        {
            return m_mean_anomaly;
        }
        else
        {
            return Util::DegreesToRadians(m_mean_anomaly);
        }
    }

    /**
     * Get the mean motion
     * @returns the mean motion (revolutions per day)
     */
    double MeanMotion() const
    {
        return m_mean_motion;
    }

    /**
     * Get the orbit number
     * @returns the orbit number
     */
    unsigned int OrbitNumber() const
    {
        return m_orbit_number;
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
        unsigned int norad_number,
        const std::string& int_designator,
        const DateTime& epoch,
        double mean_motion_dt2,
        double mean_motion_ddt6,
        double bstar,
        double inclination,
        double right_ascending_node,
        double eccentricity,
        double argument_perigee,
        double mean_anomaly,
        double mean_motion,
        unsigned int orbit_number);

    void Initialize();
    static bool IsValidLineLength(const std::string& str);
    void ExtractInteger(const std::string& str, unsigned int& val);
    void ExtractDouble(const std::string& str, int point_pos, double& val);
    void ExtractExponential(const std::string& str, double& val);

private:
    std::string m_name;
    std::string m_line_one;
    std::string m_line_two;

    std::string m_int_designator;
    DateTime m_epoch;
    double m_mean_motion_dt2{};
    double m_mean_motion_ddt6{};
    double m_bstar{};
    double m_inclination{};
    double m_right_ascending_node{};
    double m_eccentricity{};
    double m_argument_perigee{};
    double m_mean_anomaly{};
    double m_mean_motion{};
    unsigned int m_norad_number{};
    unsigned int m_orbit_number{};

    static const unsigned int TLE_LEN_LINE_DATA = 69;
    static const unsigned int TLE_LEN_LINE_NAME = 22;
};


inline std::ostream& operator<<(std::ostream& strm, const Tle& t)
{
    return strm << t.ToString();
}

} // namespace libsgp4
