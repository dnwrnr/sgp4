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


#include "OrbitalElements.h"

#include "Tle.h"

namespace libsgp4
{

OrbitalElements::OrbitalElements(const Tle& tle)
{
    /*
     * extract and format tle data
     */
    mMeanAnomaly = tle.MeanAnomaly(false);
    mAscendingNode = tle.RightAscendingNode(false);
    mArgumentPerigee = tle.ArgumentPerigee(false);
    mEccentricity = tle.Eccentricity();
    mInclination = tle.Inclination(false);
    mMeanMotion = tle.MeanMotion() * kTWOPI / kMINUTES_PER_DAY;
    mBstar = tle.BStar();
    mEpoch = tle.Epoch();

    /*
     * recover original mean motion (xnodp) and semimajor axis (aodp)
     * from input elements
     */
    const double a1 = pow(kXKE / MeanMotion(), kTWOTHIRD);
    const double cosio = cos(Inclination());
    const double theta2 = cosio * cosio;
    const double x3thm1 = 3.0 * theta2 - 1.0;
    const double eosq = Eccentricity() * Eccentricity();
    const double betao2 = 1.0 - eosq;
    const double betao = sqrt(betao2);
    const double temp = (1.5 * kCK2) * x3thm1 / (betao * betao2);
    const double del1 = temp / (a1 * a1);
    const double a0 = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + del1 * 134.0 / 81.0)));
    const double del0 = temp / (a0 * a0);

    mRecoveredMeanMotion = MeanMotion() / (1.0 + del0);
    /*
     * alternative way to calculate
     * doesnt affect final results
     * mRecoveredSemiMajorAxis = pow(XKE / RecoveredMeanMotion(), TWOTHIRD);
     */
    mRecoveredSemiMajorAxis = a0 / (1.0 - del0);

    /*
     * find perigee and period
     */
    mPerigee = (RecoveredSemiMajorAxis() * (1.0 - Eccentricity()) - kAE) * kXKMPER;
    mPeriod = kTWOPI / RecoveredMeanMotion();
}

} // namespace libsgp4
