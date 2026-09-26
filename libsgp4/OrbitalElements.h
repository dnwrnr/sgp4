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
#include "Util.h"

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
        return mMeanAnomaly;
    }

    /*
     * XNODEO
     */
    double AscendingNode() const
    {
        return mAscendingNode;
    }

    /*
     * OMEGAO
     */
    double ArgumentPerigee() const
    {
        return mArgumentPerigee;
    }

    /*
     * EO
     */
    double Eccentricity() const
    {
        return mEccentricity;
    }

    /*
     * XINCL
     */
    double Inclination() const
    {
        return mInclination;
    }

    /*
     * XNO
     */
    double MeanMotion() const
    {
        return mMeanMotion;
    }

    /*
     * BSTAR
     */
    double BStar() const
    {
        return mBstar;
    }

    /*
     * AODP
     */
    double RecoveredSemiMajorAxis() const
    {
        return mRecoveredSemiMajorAxis;
    }

    /*
     * XNODP
     */
    double RecoveredMeanMotion() const
    {
        return mRecoveredMeanMotion;
    }

    /*
     * PERIGE
     */
    double Perigee() const
    {
        return mPerigee;
    }

    /*
     * Period in minutes
     */
    double Period() const
    {
        return mPeriod;
    }

    /*
     * EPOCH
     */
    DateTime Epoch() const
    {
        return mEpoch;
    }

private:
    double mMeanAnomaly;
    double mAscendingNode;
    double mArgumentPerigee;
    double mEccentricity;
    double mInclination;
    double mMeanMotion;
    double mBstar;
    double mRecoveredSemiMajorAxis;
    double mRecoveredMeanMotion;
    double mPerigee;
    double mPeriod;
    DateTime mEpoch;
};

} // namespace libsgp4
