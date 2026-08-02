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

namespace libsgp4
{

class Tle;

/**
 * @brief The extracted orbital elements used by the SGP4 propagator.
 */
class OrbitalElements
{
public:
    explicit OrbitalElements(const Tle& tle);

    /*
     * XMO
     */
    double MeanAnomaly() const
    {
        return m_mean_anomaly;
    }

    /*
     * XNODEO
     */
    double AscendingNode() const
    {
        return m_ascending_node;
    }

    /*
     * OMEGAO
     */
    double ArgumentPerigee() const
    {
        return m_argument_perigee;
    }

    /*
     * EO
     */
    double Eccentricity() const
    {
        return m_eccentricity;
    }

    /*
     * XINCL
     */
    double Inclination() const
    {
        return m_inclination;
    }

    /*
     * XNO
     */
    double MeanMotion() const
    {
        return m_mean_motion;
    }

    /*
     * BSTAR
     */
    double BStar() const
    {
        return m_bstar;
    }

    /*
     * AODP
     */
    double RecoveredSemiMajorAxis() const
    {
        return m_recovered_semi_major_axis;
    }

    /*
     * XNODP
     */
    double RecoveredMeanMotion() const
    {
        return m_recovered_mean_motion;
    }

    /*
     * PERIGE
     */
    double Perigee() const
    {
        return m_perigee;
    }

    /*
     * Period in minutes
     */
    double Period() const
    {
        return m_period;
    }

    /*
     * EPOCH
     */
    DateTime Epoch() const
    {
        return m_epoch;
    }

private:
    double m_mean_anomaly;
    double m_ascending_node;
    double m_argument_perigee;
    double m_eccentricity;
    double m_inclination;
    double m_mean_motion;
    double m_bstar;
    double m_recovered_semi_major_axis;
    double m_recovered_mean_motion;
    double m_perigee;
    double m_period;
    DateTime m_epoch;
};

} // namespace libsgp4
