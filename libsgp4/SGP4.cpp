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

#include "SGP4.h"

#include "DecayedException.h"
#include "SatelliteException.h"
#include "Util.h"
#include "Vector.h"

#include <cmath>
#include <iomanip>

namespace libsgp4
{

void SGP4::SetTle(const Tle& tle)
{
    /*
     * extract and format tle data
     */
    mElements = OrbitalElements(tle);

    Initialise();
}

/*
 * traceability: short identifier -> descriptive name
  aycof                    -> ayCoeff                            betao                    -> sqrtOneMinusEccSq
  betao2                   -> oneMinusEccSq                      c1                       -> dragCoeff
  c1sq                     -> dragCoeffSquared                   c2                       -> dragCoeff2
  c3                       -> dragCoeff3                         c4                       -> dragCoeff4
  c5                       -> dragCoeff5                         coef                     -> atmCoef
  coef1                    -> atmCoef1                           cosio                    -> cosInclination0
  d2                       -> d2Coeff                            d3                       -> d3Coeff
  d4                       -> d4Coeff                            delmo                    -> deltaMeanAnomaly0
  eeta                     -> eccTimesEta                        eosq                     -> eccentricitySq
  etasq                    -> etaSquared                         gsto                     -> greenwichSiderealTime
  omgcof                   -> argPerigeeDragCoeff                omgdot                   -> argPerigeeDot
  pinvsq                   -> pinvSquared                        psisq                    -> psiSquared
  qoms24                   -> qoms2t                             s4                       -> densityHeightS
  sinio                    -> sinInclination0                    sinmo                    -> sinMeanAnomaly0
  t2cof                    -> t2Coeff                            t3cof                    -> t3Coeff
  t4cof                    -> t4Coeff                            t5cof                    -> t5Coeff
  temp                     -> d3PreFactor                        temp1                    -> ck2Term1
  temp2                    -> ck2Term2                           temp3                    -> ck4Term
  theta2                   -> cosSqInc                           theta4                   -> cos4Inc
  tsi                      -> inverseSmaMinusS                   x1m5th                   -> oneMinus5CosSqInc
  x1mth2                   -> sinSqInc                           x3thm1                   -> threeCosSqIncMinus1
  x7thm1                   -> sevenCosSqIncMinus1                xhdot1                   -> raanDotFirstOrder
  xlcof                    -> xlCoeff                            xmcof                    -> meanAnomalyDragCoeff
  xmdot                    -> meanAnomalyDot                     xnodcf                   -> raanDragCoeff
  xnodot                   -> raanDot
 */
void SGP4::Initialise()
{
    /*
     * reset all constants etc
     */
    Reset();

    /*
     * error checks
     */
    if (mElements.Eccentricity() < 0.0 || mElements.Eccentricity() > 0.999)
    {
        throw SatelliteException("Eccentricity out of range");
    }

    if (mElements.Inclination() < 0.0 || mElements.Inclination() > kPI)
    {
        throw SatelliteException("Inclination out of range");
    }

    RecomputeConstants(mElements.Inclination(),
                       mCommonConsts.sinInclination0,
                       mCommonConsts.cosInclination0,
                       mCommonConsts.threeCosSqIncMinus1,
                       mCommonConsts.sinSqInc,
                       mCommonConsts.sevenCosSqIncMinus1,
                       mCommonConsts.xlCoeff,
                       mCommonConsts.ayCoeff);

    const double cosSqInc = mCommonConsts.cosInclination0 * mCommonConsts.cosInclination0;
    const double eccentricitySq = mElements.Eccentricity() * mElements.Eccentricity();
    const double oneMinusEccSq = 1.0 - eccentricitySq;
    const double sqrtOneMinusEccSq = sqrt(oneMinusEccSq);

    if (mElements.Period() >= 225.0)
    {
        mUseDeepSpace = true;
    }
    else
    {
        mUseDeepSpace = false;
        mUseSimpleModel = false;
        /*
         * for perigee less than 220 kilometers, the simpleModel flag is set
         * and the equations are truncated to linear variation in sqrt a and
         * quadratic variation in mean anomly. also, the c3 term, the
         * delta omega term and the delta m term are dropped
         */
        if (mElements.Perigee() < 220.0)
        {
            mUseSimpleModel = true;
        }
    }

    /*
     * for perigee below 156km, the values of
     * s4 and qoms2t are altered
     */
    double densityHeightS = kS;
    double qoms2t = kQOMS2T;
    if (mElements.Perigee() < 156.0)
    {
        densityHeightS = mElements.Perigee() - 78.0;
        if (mElements.Perigee() < 98.0)
        {
            densityHeightS = 20.0;
        }
        qoms2t = pow((120.0 - densityHeightS) * kAE / kXKMPER, 4.0);
        densityHeightS = densityHeightS / kXKMPER + kAE;
    }

    /*
     * generate constants
     */
    const double pinvSquared =
        1.0 / (mElements.RecoveredSemiMajorAxis() * mElements.RecoveredSemiMajorAxis() * oneMinusEccSq * oneMinusEccSq);
    const double inverseSmaMinusS = 1.0 / (mElements.RecoveredSemiMajorAxis() - densityHeightS);
    mCommonConsts.eta = mElements.RecoveredSemiMajorAxis() * mElements.Eccentricity() * inverseSmaMinusS;
    const double etaSquared = mCommonConsts.eta * mCommonConsts.eta;
    const double eccTimesEta = mElements.Eccentricity() * mCommonConsts.eta;
    const double psiSquared = std::abs(1.0 - etaSquared);
    const double atmCoef = qoms2t * pow(inverseSmaMinusS, 4.0);
    const double atmCoef1 = atmCoef / pow(psiSquared, 3.5);
    const double dragCoeff2 =
        atmCoef1 * mElements.RecoveredMeanMotion() *
        (mElements.RecoveredSemiMajorAxis() * (1.0 + 1.5 * etaSquared + eccTimesEta * (4.0 + etaSquared)) +
         0.75 * kCK2 * inverseSmaMinusS / psiSquared * mCommonConsts.threeCosSqIncMinus1 *
             (8.0 + 3.0 * etaSquared * (8.0 + etaSquared)));
    mCommonConsts.dragCoeff = mElements.BStar() * dragCoeff2;
    mCommonConsts.dragCoeff4 =
        2.0 * mElements.RecoveredMeanMotion() * atmCoef1 * mElements.RecoveredSemiMajorAxis() * oneMinusEccSq *
        (mCommonConsts.eta * (2.0 + 0.5 * etaSquared) + mElements.Eccentricity() * (0.5 + 2.0 * etaSquared) -
         2.0 * kCK2 * inverseSmaMinusS / (mElements.RecoveredSemiMajorAxis() * psiSquared) *
             (-3.0 * mCommonConsts.threeCosSqIncMinus1 *
                  (1.0 - 2.0 * eccTimesEta + etaSquared * (1.5 - 0.5 * eccTimesEta)) +
              0.75 * mCommonConsts.sinSqInc * (2.0 * etaSquared - eccTimesEta * (1.0 + etaSquared)) *
                  cos(2.0 * mElements.ArgumentPerigee())));
    const double cos4Inc = cosSqInc * cosSqInc;
    const double ck2Term1 = 3.0 * kCK2 * pinvSquared * mElements.RecoveredMeanMotion();
    const double ck2Term2 = ck2Term1 * kCK2 * pinvSquared;
    const double ck4Term = 1.25 * kCK4 * pinvSquared * pinvSquared * mElements.RecoveredMeanMotion();
    mCommonConsts.meanAnomalyDot = mElements.RecoveredMeanMotion() +
                                   0.5 * ck2Term1 * sqrtOneMinusEccSq * mCommonConsts.threeCosSqIncMinus1 +
                                   0.0625 * ck2Term2 * sqrtOneMinusEccSq * (13.0 - 78.0 * cosSqInc + 137.0 * cos4Inc);
    const double oneMinus5CosSqInc = 1.0 - 5.0 * cosSqInc;
    mCommonConsts.argPerigeeDot = -0.5 * ck2Term1 * oneMinus5CosSqInc +
                                  0.0625 * ck2Term2 * (7.0 - 114.0 * cosSqInc + 395.0 * cos4Inc) +
                                  ck4Term * (3.0 - 36.0 * cosSqInc + 49.0 * cos4Inc);
    const double raanDotFirstOrder = -ck2Term1 * mCommonConsts.cosInclination0;
    mCommonConsts.raanDot =
        raanDotFirstOrder + (0.5 * ck2Term2 * (4.0 - 19.0 * cosSqInc) + 2.0 * ck4Term * (3.0 - 7.0 * cosSqInc)) *
                                mCommonConsts.cosInclination0;
    mCommonConsts.raanDragCoeff = 3.5 * oneMinusEccSq * raanDotFirstOrder * mCommonConsts.dragCoeff;
    mCommonConsts.t2Coeff = 1.5 * mCommonConsts.dragCoeff;

    if (mUseDeepSpace)
    {
        mDeepspaceConsts.greenwichSiderealTime = mElements.Epoch().ToGreenwichSiderealTime();

        DeepSpaceInitialise(eccentricitySq,
                            mCommonConsts.sinInclination0,
                            mCommonConsts.cosInclination0,
                            sqrtOneMinusEccSq,
                            cosSqInc,
                            oneMinusEccSq,
                            mCommonConsts.meanAnomalyDot,
                            mCommonConsts.argPerigeeDot,
                            mCommonConsts.raanDot);
    }
    else
    {
        double dragCoeff3 = 0.0;
        if (mElements.Eccentricity() > 1.0e-4)
        {
            dragCoeff3 = atmCoef * inverseSmaMinusS * kA3OVK2 * mElements.RecoveredMeanMotion() * kAE *
                         mCommonConsts.sinInclination0 / mElements.Eccentricity();
        }

        mNearspaceConsts.dragCoeff5 = 2.0 * atmCoef1 * mElements.RecoveredSemiMajorAxis() * oneMinusEccSq *
                                      (1.0 + 2.75 * (etaSquared + eccTimesEta) + eccTimesEta * etaSquared);
        mNearspaceConsts.argPerigeeDragCoeff = mElements.BStar() * dragCoeff3 * cos(mElements.ArgumentPerigee());

        mNearspaceConsts.meanAnomalyDragCoeff = 0.0;
        if (mElements.Eccentricity() > 1.0e-4)
        {
            mNearspaceConsts.meanAnomalyDragCoeff = -kTWOTHIRD * atmCoef * mElements.BStar() * kAE / eccTimesEta;
        }

        mNearspaceConsts.deltaMeanAnomaly0 = pow(1.0 + mCommonConsts.eta * (cos(mElements.MeanAnomaly())), 3.0);
        mNearspaceConsts.sinMeanAnomaly0 = sin(mElements.MeanAnomaly());

        if (!mUseSimpleModel)
        {
            const double dragCoeffSquared = mCommonConsts.dragCoeff * mCommonConsts.dragCoeff;
            mNearspaceConsts.d2Coeff = 4.0 * mElements.RecoveredSemiMajorAxis() * inverseSmaMinusS * dragCoeffSquared;
            const double d3PreFactor = mNearspaceConsts.d2Coeff * inverseSmaMinusS * mCommonConsts.dragCoeff / 3.0;
            mNearspaceConsts.d3Coeff = (17.0 * mElements.RecoveredSemiMajorAxis() + densityHeightS) * d3PreFactor;
            mNearspaceConsts.d4Coeff = 0.5 * d3PreFactor * mElements.RecoveredSemiMajorAxis() * inverseSmaMinusS *
                                       (221.0 * mElements.RecoveredSemiMajorAxis() + 31.0 * densityHeightS) *
                                       mCommonConsts.dragCoeff;
            mNearspaceConsts.t3Coeff = mNearspaceConsts.d2Coeff + 2.0 * dragCoeffSquared;
            mNearspaceConsts.t4Coeff =
                0.25 * (3.0 * mNearspaceConsts.d3Coeff +
                        mCommonConsts.dragCoeff * (12.0 * mNearspaceConsts.d2Coeff + 10.0 * dragCoeffSquared));
            mNearspaceConsts.t5Coeff =
                0.2 * (3.0 * mNearspaceConsts.d4Coeff + 12.0 * mCommonConsts.dragCoeff * mNearspaceConsts.d3Coeff +
                       6.0 * mNearspaceConsts.d2Coeff * mNearspaceConsts.d2Coeff +
                       15.0 * dragCoeffSquared * (2.0 * mNearspaceConsts.d2Coeff + dragCoeffSquared));
        }
    }
}

Eci SGP4::FindPosition(const DateTime& dt) const
{
    return FindPosition((dt - mElements.Epoch()).TotalMinutes());
}

Eci SGP4::FindPosition(double tsince) const
{
    if (!std::isfinite(tsince))
    {
        throw SatelliteException("Error: (tsince not finite)");
    }

    if (mUseDeepSpace)
    {
        return FindPositionSDP4(tsince);
    }
    else
    {
        return FindPositionSGP4(tsince);
    }
}

/*
 * traceability: short identifier -> descriptive name
  a                        -> semiMajorAxis                      c1                       -> dragCoeff
  c4                       -> dragCoeff4                         e                        -> eccentricity
  em                       -> eccentricityFromEpoch              omega                    -> argPerigee
  omgadf                   -> argPerigeeSecular                  omgdot                   -> argPerigeeDot
  perturbedAycof           -> perturbedAyCoeff                   perturbedCosio           -> perturbedCosInclination0
  perturbedSinio           -> perturbedSinInclination0           perturbedX1mth2          -> perturbedSinSqInc
  perturbedX3thm1          -> perturbedThreeCosSqIncMinus1       perturbedX7thm1          ->
 perturbedSevenCosSqIncMinus1 perturbedXlcof           -> perturbedXlCoeff                   t2cof                    ->
 t2Coeff tempa                    -> semiMajorAxisDrag                  tempe                    -> eccentricityDrag
  templ                    -> meanLongitudeDrag                  tsq                      -> tsinceSquared
  v                        -> smaFactor                          xinc                     -> inclination
  xl                       -> meanLongitude                      xmam                     -> meanAnomaly
  xmdf                     -> meanAnomalySecular                 xmdot                    -> meanAnomalyDot
  xn                       -> meanMotion                         xnodcf                   -> raanDragCoeff
  xnoddf                   -> raanSecular                        xnode                    -> raan
  xnodot                   -> raanDot
 */
Eci SGP4::FindPositionSDP4(double tsince) const
{
    /*
     * the final values
     */
    double eccentricity;
    double semiMajorAxis;
    double argPerigee;
    double meanLongitude;
    double raan;
    double inclination;

    /*
     * update for secular gravity and atmospheric drag
     */
    double meanAnomalySecular = mElements.MeanAnomaly() + mCommonConsts.meanAnomalyDot * tsince;
    double argPerigeeSecular = mElements.ArgumentPerigee() + mCommonConsts.argPerigeeDot * tsince;
    const double raanSecular = mElements.AscendingNode() + mCommonConsts.raanDot * tsince;

    const double tsinceSquared = tsince * tsince;
    raan = raanSecular + mCommonConsts.raanDragCoeff * tsinceSquared;
    double semiMajorAxisDrag = 1.0 - mCommonConsts.dragCoeff * tsince;
    double eccentricityDrag = mElements.BStar() * mCommonConsts.dragCoeff4 * tsince;
    double meanLongitudeDrag = mCommonConsts.t2Coeff * tsinceSquared;

    double meanMotion = mElements.RecoveredMeanMotion();
    double eccentricityFromEpoch = mElements.Eccentricity();
    inclination = mElements.Inclination();

    DeepSpaceSecular(tsince,
                     mElements,
                     mCommonConsts,
                     mDeepspaceConsts,
                     mIntegratorParams,
                     meanAnomalySecular,
                     argPerigeeSecular,
                     raan,
                     eccentricityFromEpoch,
                     inclination,
                     meanMotion);

    if (!(meanMotion > 0.0))
    {
        throw SatelliteException("Error: (xn <= 0.0 or not finite)");
    }

    const double smaFactor = kXKE / meanMotion;
    semiMajorAxis = cbrt(smaFactor * smaFactor) * semiMajorAxisDrag * semiMajorAxisDrag;
    eccentricity = eccentricityFromEpoch - eccentricityDrag;
    double meanAnomaly = meanAnomalySecular + mElements.RecoveredMeanMotion() * meanLongitudeDrag;

    DeepSpacePeriodics(tsince, mDeepspaceConsts, eccentricity, inclination, argPerigeeSecular, raan, meanAnomaly);

    /*
     * keeping xinc positive important unless you need to display xinc
     * and dislike negative inclinations
     */
    if (inclination < 0.0)
    {
        inclination = -inclination;
        raan += kPI;
        argPerigeeSecular -= kPI;
    }

    meanLongitude = meanAnomaly + argPerigeeSecular + raan;
    argPerigee = argPerigeeSecular;

    /*
     * fix tolerance for error recognition
     */
    if (!(eccentricity > -0.001))
    {
        throw SatelliteException("Error: (e <= -0.001 or not finite)");
    }
    else if (eccentricity < 1.0e-6)
    {
        eccentricity = 1.0e-6;
    }
    else if (eccentricity > (1.0 - 1.0e-6))
    {
        eccentricity = 1.0 - 1.0e-6;
    }

    /*
     * re-compute the perturbed values
     */
    double perturbedSinInclination0;
    double perturbedCosInclination0;
    double perturbedThreeCosSqIncMinus1;
    double perturbedSinSqInc;
    double perturbedSevenCosSqIncMinus1;
    double perturbedXlCoeff;
    double perturbedAyCoeff;
    RecomputeConstants(inclination,
                       perturbedSinInclination0,
                       perturbedCosInclination0,
                       perturbedThreeCosSqIncMinus1,
                       perturbedSinSqInc,
                       perturbedSevenCosSqIncMinus1,
                       perturbedXlCoeff,
                       perturbedAyCoeff);

    /*
     * using calculated values, find position and velocity
     */
    return CalculateFinalPositionVelocity(mElements.Epoch().AddMinutes(tsince),
                                          eccentricity,
                                          semiMajorAxis,
                                          argPerigee,
                                          meanLongitude,
                                          raan,
                                          inclination,
                                          perturbedXlCoeff,
                                          perturbedAyCoeff,
                                          perturbedThreeCosSqIncMinus1,
                                          perturbedSinSqInc,
                                          perturbedSevenCosSqIncMinus1,
                                          perturbedCosInclination0,
                                          perturbedSinInclination0);
}

/*
 * traceability: short identifier -> descriptive name
  aycof                    -> ayCoeff                            cosio                    -> cosInclination0
  sinio                    -> sinInclination0                    theta2                   -> cosSqInc
  x1mth2                   -> sinSqInc                           x3thm1                   -> threeCosSqIncMinus1
  x7thm1                   -> sevenCosSqIncMinus1                xinc                     -> inclination
  xlcof                    -> xlCoeff
 */
void SGP4::RecomputeConstants(double inclination,
                              double& sinInclination0,
                              double& cosInclination0,
                              double& threeCosSqIncMinus1,
                              double& sinSqInc,
                              double& sevenCosSqIncMinus1,
                              double& xlCoeff,
                              double& ayCoeff)
{
    sinInclination0 = sin(inclination);
    cosInclination0 = cos(inclination);

    const double cosSqInc = cosInclination0 * cosInclination0;

    threeCosSqIncMinus1 = 3.0 * cosSqInc - 1.0;
    sinSqInc = 1.0 - cosSqInc;
    sevenCosSqIncMinus1 = 7.0 * cosSqInc - 1.0;

    if (std::abs(cosInclination0 + 1.0) > 1.5e-12)
    {
        xlCoeff = 0.125 * kA3OVK2 * sinInclination0 * (3.0 + 5.0 * cosInclination0) / (1.0 + cosInclination0);
    }
    else
    {
        xlCoeff = 0.125 * kA3OVK2 * sinInclination0 * (3.0 + 5.0 * cosInclination0) / 1.5e-12;
    }

    ayCoeff = 0.25 * kA3OVK2 * sinInclination0;
}

/*
 * traceability: short identifier -> descriptive name
  a                        -> semiMajorAxis                      aycof                    -> ayCoeff
  c1                       -> dragCoeff                          c4                       -> dragCoeff4
  c5                       -> dragCoeff5                         cosio                    -> cosInclination0
  d2                       -> d2Coeff                            d3                       -> d3Coeff
  d4                       -> d4Coeff                            delm                     -> meanAnomalyCorrection
  delmo                    -> deltaMeanAnomaly0                  delomg                   -> argPerigeeCorrection
  e                        -> eccentricity                       omega                    -> argPerigee
  omgadf                   -> argPerigeeSecular                  omgcof                   -> argPerigeeDragCoeff
  omgdot                   -> argPerigeeDot                      sinio                    -> sinInclination0
  sinmo                    -> sinMeanAnomaly0                    t2cof                    -> t2Coeff
  t3cof                    -> t3Coeff                            t4cof                    -> t4Coeff
  t5cof                    -> t5Coeff                            tcube                    -> tsinceCubed
  temp                     -> combinedCorrection                 tempa                    -> semiMajorAxisDrag
  tempe                    -> eccentricityDrag                   templ                    -> meanLongitudeDrag
  tfour                    -> tsinceFourth                       tsq                      -> tsinceSquared
  x1mth2                   -> sinSqInc                           x1p                      -> onePlusEtaCosM
  x3thm1                   -> threeCosSqIncMinus1                x7thm1                   -> sevenCosSqIncMinus1
  xinc                     -> inclination                        xl                       -> meanLongitude
  xlcof                    -> xlCoeff                            xmcof                    -> meanAnomalyDragCoeff
  xmdf                     -> meanAnomalySecular                 xmdot                    -> meanAnomalyDot
  xmp                      -> meanAnomalyPerturbed               xnodcf                   -> raanDragCoeff
  xnoddf                   -> raanSecular                        xnode                    -> raan
  xnodot                   -> raanDot
 */
Eci SGP4::FindPositionSGP4(double tsince) const
{
    /*
     * the final values
     */
    double eccentricity;
    double semiMajorAxis;
    double argPerigee;
    double meanLongitude;
    double raan;
    const double inclination = mElements.Inclination();

    /*
     * update for secular gravity and atmospheric drag
     */
    const double meanAnomalySecular = mElements.MeanAnomaly() + mCommonConsts.meanAnomalyDot * tsince;
    const double argPerigeeSecular = mElements.ArgumentPerigee() + mCommonConsts.argPerigeeDot * tsince;
    const double raanSecular = mElements.AscendingNode() + mCommonConsts.raanDot * tsince;

    argPerigee = argPerigeeSecular;
    double meanAnomalyPerturbed = meanAnomalySecular;

    const double tsinceSquared = tsince * tsince;
    raan = raanSecular + mCommonConsts.raanDragCoeff * tsinceSquared;
    double semiMajorAxisDrag = 1.0 - mCommonConsts.dragCoeff * tsince;
    double eccentricityDrag = mElements.BStar() * mCommonConsts.dragCoeff4 * tsince;
    double meanLongitudeDrag = mCommonConsts.t2Coeff * tsinceSquared;

    if (!mUseSimpleModel)
    {
        const double argPerigeeCorrection = mNearspaceConsts.argPerigeeDragCoeff * tsince;
        const double onePlusEtaCosM = 1.0 + mCommonConsts.eta * cos(Util::WrapTwoPI(meanAnomalySecular));
        const double meanAnomalyCorrection =
            mNearspaceConsts.meanAnomalyDragCoeff *
            (onePlusEtaCosM * onePlusEtaCosM * onePlusEtaCosM - mNearspaceConsts.deltaMeanAnomaly0);
        const double combinedCorrection = argPerigeeCorrection + meanAnomalyCorrection;

        meanAnomalyPerturbed += combinedCorrection;
        argPerigee -= combinedCorrection;

        const double tsinceCubed = tsinceSquared * tsince;
        const double tsinceFourth = tsince * tsinceCubed;

        semiMajorAxisDrag = semiMajorAxisDrag - mNearspaceConsts.d2Coeff * tsinceSquared -
                            mNearspaceConsts.d3Coeff * tsinceCubed - mNearspaceConsts.d4Coeff * tsinceFourth;
        eccentricityDrag += mElements.BStar() * mNearspaceConsts.dragCoeff5 *
                            (sin(Util::WrapTwoPI(meanAnomalyPerturbed)) - mNearspaceConsts.sinMeanAnomaly0);
        meanLongitudeDrag += mNearspaceConsts.t3Coeff * tsinceCubed +
                             tsinceFourth * (mNearspaceConsts.t4Coeff + tsince * mNearspaceConsts.t5Coeff);
    }

    semiMajorAxis = mElements.RecoveredSemiMajorAxis() * semiMajorAxisDrag * semiMajorAxisDrag;
    eccentricity = mElements.Eccentricity() - eccentricityDrag;
    meanLongitude = meanAnomalyPerturbed + argPerigee + raan + mElements.RecoveredMeanMotion() * meanLongitudeDrag;

    /*
     * fix tolerance for error recognition
     */
    if (!(eccentricity > -0.001))
    {
        throw SatelliteException("Error: (e <= -0.001 or not finite)");
    }
    else if (eccentricity < 1.0e-6)
    {
        eccentricity = 1.0e-6;
    }
    else if (eccentricity > (1.0 - 1.0e-6))
    {
        eccentricity = 1.0 - 1.0e-6;
    }

    /*
     * using calculated values, find position and velocity
     * we can pass in constants from Initialise() as these dont change
     */
    return CalculateFinalPositionVelocity(mElements.Epoch().AddMinutes(tsince),
                                          eccentricity,
                                          semiMajorAxis,
                                          argPerigee,
                                          meanLongitude,
                                          raan,
                                          inclination,
                                          mCommonConsts.xlCoeff,
                                          mCommonConsts.ayCoeff,
                                          mCommonConsts.threeCosSqIncMinus1,
                                          mCommonConsts.sinSqInc,
                                          mCommonConsts.sevenCosSqIncMinus1,
                                          mCommonConsts.cosInclination0,
                                          mCommonConsts.sinInclination0);
}

/*
 * traceability: short identifier -> descriptive name
  a                        -> semiMajorAxis                      axn                      -> eCosArgPerigee
  aycof                    -> ayCoeff                            ayn                      -> eSinArgPerigee
  aynl                     -> longPeriodAy                       beta2                    -> oneMinusEccSq
  betal                    -> sqrtOneMinusEVectorSq              capu                     -> keplerMeanAnomaly
  cos2u                    -> cosTwoArgLat                       cosdelta                 -> cosDelta
  cosepw                   -> cosEccentricAnomaly                cosik                    -> cosIk
  cosio                    -> cosInclination0                    cosnok                   -> cosNodeK
  cosu                     -> cosArgLatitude                     cosuk                    -> cosUk
  delta                    -> argLatitudePerturbation            deltaEpw                 -> eccentricAnomalyCorrection
  e                        -> eccentricity                       ecose                    -> eCosE
  elsq                     -> eVectorSq                          epw                      -> eccentricAnomaly
  esine                    -> eSinE                              f                        -> keplerResidual
  fdot                     -> keplerResidualDerivative           keplerRunning            -> keplerConverged
  maxNewtonNaphson         -> keplerStepLimit                    omega                    -> argPerigee
  pl                       -> semiLatusRectum                    r                        -> radius
  rdot                     -> radialVelocity                     rdotk                    -> radialVelocityPerturbed
  rfdot                    -> transverseVelocity                 rfdotk                   -> transverseVelocityPerturbed
  rk                       -> radiusPerturbed                    sin2u                    -> sinTwoArgLat
  sindelta                 -> sinDelta                           sinepw                   -> sinEccentricAnomaly
  sinik                    -> sinIk                              sinio                    -> sinInclination0
  sinnok                   -> sinNodeK                           sinu                     -> sinArgLatitude
  sinuk                    -> sinUk                              sqrtA                    -> sqrtSemiMajorAxis
  temp11                   -> oneOverSmaOneMinusEccSq            temp21                   -> oneMinusEVectorSq
  temp31                   -> oneOverRadius                      temp32                   -> smaOverRadius
  temp33                   -> oneOverOnePlusBeta                 temp41                   -> oneOverSemiLatusRectum
  temp42                   -> ck2OverP                           temp43                   -> ck2OverPSq
  ux                       -> uBasisX                            uy                       -> uBasisY
  uz                       -> uBasisZ                            vx                       -> vBasisX
  vy                       -> vBasisY                            vz                       -> vBasisZ
  x                        -> posX                               x1mth2                   -> sinSqInc
  x3thm1                   -> threeCosSqIncMinus1                x7thm1                   -> sevenCosSqIncMinus1
  xdot                     -> velX                               xinc                     -> inclination
  xinck                    -> inclinationPerturbed               xl                       -> meanLongitude
  xlcof                    -> xlCoeff                            xll                      -> longPeriodCorrection
  xlt                      -> longitudePerturbed                 xmx                      -> nodeVectorX
  xmy                      -> nodeVectorY                        xn                       -> meanMotion
  xnode                    -> raan                               xnodek                   -> raanPerturbed
  y                        -> posY                               ydot                     -> velY
  z                        -> posZ                               zdot                     -> velZ
 */
Eci SGP4::CalculateFinalPositionVelocity(const DateTime& dt,
                                         double eccentricity,
                                         double semiMajorAxis,
                                         double argPerigee,
                                         double meanLongitude,
                                         double raan,
                                         double inclination,
                                         double xlCoeff,
                                         double ayCoeff,
                                         double threeCosSqIncMinus1,
                                         double sinSqInc,
                                         double sevenCosSqIncMinus1,
                                         double cosInclination0,
                                         double sinInclination0)
{
    if (!std::isfinite(eccentricity) || !std::isfinite(semiMajorAxis) || !std::isfinite(argPerigee) ||
        !std::isfinite(meanLongitude) || !std::isfinite(raan) || !std::isfinite(inclination))
    {
        throw SatelliteException("Error: (element not finite)");
    }

    const double oneMinusEccSq = 1.0 - eccentricity * eccentricity;
    const double sqrtSemiMajorAxis = sqrt(semiMajorAxis);
    const double meanMotion = kXKE / (semiMajorAxis * sqrtSemiMajorAxis);
    /*
     * long period periodics
     */
    const double eCosArgPerigee = eccentricity * cos(argPerigee);
    const double oneOverSmaOneMinusEccSq = 1.0 / (semiMajorAxis * oneMinusEccSq);
    const double longPeriodCorrection = oneOverSmaOneMinusEccSq * xlCoeff * eCosArgPerigee;
    const double longPeriodAy = oneOverSmaOneMinusEccSq * ayCoeff;
    const double longitudePerturbed = meanLongitude + longPeriodCorrection;
    const double eSinArgPerigee = eccentricity * sin(argPerigee) + longPeriodAy;
    const double eVectorSq = eCosArgPerigee * eCosArgPerigee + eSinArgPerigee * eSinArgPerigee;

    if (!(eVectorSq < 1.0))
    {
        throw SatelliteException("Error: (elsq >= 1.0 or not finite)");
    }

    /*
     * solve keplers equation
     * - solve using Newton-Raphson root solving
     * - here capu is almost the mean anomoly
     * - initialise the eccentric anomaly term epw
     * - The fmod saves reduction of angle to +/-2pi in sin/cos() and prevents
     * convergence problems.
     */
    const double keplerMeanAnomaly = fmod(longitudePerturbed - raan, kTWOPI);
    double eccentricAnomaly = keplerMeanAnomaly;

    double sinEccentricAnomaly = 0.0;
    double cosEccentricAnomaly = 0.0;
    double eCosE = 0.0;
    double eSinE = 0.0;

    /*
     * sensibility check for N-R correction
     */
    const double keplerStepLimit = 1.25 * std::abs(sqrt(eVectorSq));

    bool keplerConverged = true;

    for (int i = 0; i < 10 && keplerConverged; i++)
    {
        sinEccentricAnomaly = sin(eccentricAnomaly);
        cosEccentricAnomaly = cos(eccentricAnomaly);
        eCosE = eCosArgPerigee * cosEccentricAnomaly + eSinArgPerigee * sinEccentricAnomaly;
        eSinE = eCosArgPerigee * sinEccentricAnomaly - eSinArgPerigee * cosEccentricAnomaly;

        double keplerResidual = keplerMeanAnomaly - eccentricAnomaly + eSinE;

        if (std::abs(keplerResidual) < 1.0e-12)
        {
            keplerConverged = false;
        }
        else
        {
            /*
             * 1st order Newton-Raphson correction
             */
            const double keplerResidualDerivative = 1.0 - eCosE;
            double eccentricAnomalyCorrection = keplerResidual / keplerResidualDerivative;

            /*
             * 2nd order Newton-Raphson correction.
             * f / (fdot - 0.5 * d2f * f/fdot)
             */
            if (i == 0)
            {
                if (eccentricAnomalyCorrection > keplerStepLimit)
                {
                    eccentricAnomalyCorrection = keplerStepLimit;
                }
                else if (eccentricAnomalyCorrection < -keplerStepLimit)
                {
                    eccentricAnomalyCorrection = -keplerStepLimit;
                }
            }
            else
            {
                eccentricAnomalyCorrection =
                    keplerResidual / (keplerResidualDerivative + 0.5 * eSinE * eccentricAnomalyCorrection);
            }

            /*
             * Newton-Raphson correction of -F/DF
             */
            eccentricAnomaly += eccentricAnomalyCorrection;
        }
    }
    /*
     * short period preliminary quantities
     */
    const double oneMinusEVectorSq = 1.0 - eVectorSq;
    const double semiLatusRectum = semiMajorAxis * oneMinusEVectorSq;

    if (!(semiLatusRectum >= 0.0))
    {
        throw SatelliteException("Error: (pl < 0.0 or not finite)");
    }

    const double radius = semiMajorAxis * (1.0 - eCosE);
    const double oneOverRadius = 1.0 / radius;
    const double radialVelocity = kXKE * sqrtSemiMajorAxis * eSinE * oneOverRadius;
    const double transverseVelocity = kXKE * sqrt(semiLatusRectum) * oneOverRadius;
    const double smaOverRadius = semiMajorAxis * oneOverRadius;
    const double sqrtOneMinusEVectorSq = sqrt(oneMinusEVectorSq);
    const double oneOverOnePlusBeta = 1.0 / (1.0 + sqrtOneMinusEVectorSq);
    const double cosArgLatitude =
        smaOverRadius * (cosEccentricAnomaly - eCosArgPerigee + eSinArgPerigee * eSinE * oneOverOnePlusBeta);
    const double sinArgLatitude =
        smaOverRadius * (sinEccentricAnomaly - eSinArgPerigee - eCosArgPerigee * eSinE * oneOverOnePlusBeta);
    const double sinTwoArgLat = 2.0 * sinArgLatitude * cosArgLatitude;
    const double cosTwoArgLat = 2.0 * cosArgLatitude * cosArgLatitude - 1.0;

    /*
     * update for short periodics
     */
    const double oneOverSemiLatusRectum = 1.0 / semiLatusRectum;
    const double ck2OverP = kCK2 * oneOverSemiLatusRectum;
    const double ck2OverPSq = ck2OverP * oneOverSemiLatusRectum;

    const double radiusPerturbed = radius * (1.0 - 1.5 * ck2OverPSq * sqrtOneMinusEVectorSq * threeCosSqIncMinus1) +
                                   0.5 * ck2OverP * sinSqInc * cosTwoArgLat;
    const double argLatitudePerturbation = -0.25 * ck2OverPSq * sevenCosSqIncMinus1 * sinTwoArgLat;
    const double raanPerturbed = raan + 1.5 * ck2OverPSq * cosInclination0 * sinTwoArgLat;
    const double inclinationPerturbed =
        inclination + 1.5 * ck2OverPSq * cosInclination0 * sinInclination0 * cosTwoArgLat;
    const double radialVelocityPerturbed = radialVelocity - meanMotion * ck2OverP * sinSqInc * sinTwoArgLat;
    const double transverseVelocityPerturbed =
        transverseVelocity + meanMotion * ck2OverP * (sinSqInc * cosTwoArgLat + 1.5 * threeCosSqIncMinus1);

    /*
     * orientation vectors
     */
    const double sinDelta = sin(argLatitudePerturbation);
    const double cosDelta = cos(argLatitudePerturbation);
    const double sinUk = sinArgLatitude * cosDelta + cosArgLatitude * sinDelta;
    const double cosUk = cosArgLatitude * cosDelta - sinArgLatitude * sinDelta;
    const double sinIk = sin(inclinationPerturbed);
    const double cosIk = cos(inclinationPerturbed);
    const double sinNodeK = sin(raanPerturbed);
    const double cosNodeK = cos(raanPerturbed);
    const double nodeVectorX = -sinNodeK * cosIk;
    const double nodeVectorY = cosNodeK * cosIk;
    const double uBasisX = nodeVectorX * sinUk + cosNodeK * cosUk;
    const double uBasisY = nodeVectorY * sinUk + sinNodeK * cosUk;
    const double uBasisZ = sinIk * sinUk;
    const double vBasisX = nodeVectorX * cosUk - cosNodeK * sinUk;
    const double vBasisY = nodeVectorY * cosUk - sinNodeK * sinUk;
    const double vBasisZ = sinIk * cosUk;
    /*
     * position and velocity
     */
    const double posX = radiusPerturbed * uBasisX * kXKMPER;
    const double posY = radiusPerturbed * uBasisY * kXKMPER;
    const double posZ = radiusPerturbed * uBasisZ * kXKMPER;
    Vector position(posX, posY, posZ);
    const double velX = (radialVelocityPerturbed * uBasisX + transverseVelocityPerturbed * vBasisX) * kXKMPER / 60.0;
    const double velY = (radialVelocityPerturbed * uBasisY + transverseVelocityPerturbed * vBasisY) * kXKMPER / 60.0;
    const double velZ = (radialVelocityPerturbed * uBasisZ + transverseVelocityPerturbed * vBasisZ) * kXKMPER / 60.0;
    Vector velocity(velX, velY, velZ);

    if (!(radiusPerturbed >= 1.0))
    {
        if (!std::isfinite(radiusPerturbed))
        {
            throw SatelliteException("Error: (position radius not finite)");
        }
        throw DecayedException(dt, position, velocity);
    }

    return Eci(dt, position, velocity);
}

static inline double EvaluateCubicPolynomial(const double x,
                                             const double constant,
                                             const double linear,
                                             const double squared,
                                             const double cubed)
{
    return constant + x * linear + x * x * squared + x * x * x * cubed;
}

/*
 * traceability: short identifier -> descriptive name
  C1L                      -> kLunarC1                           C1SS                     -> kSolarC1
  Q22                      -> kSynchronousQ22                    Q31                      -> kSynchronousQ31
  Q33                      -> kSynchronousQ33                    ROOT22                   -> kResonanceRoot22
  ROOT32                   -> kResonanceRoot32                   ROOT44                   -> kResonanceRoot44
  ROOT52                   -> kResonanceRoot52                   ROOT54                   -> kResonanceRoot54
  ZCOSGS                   -> kCosSolarArgPerigee                ZCOSIS                   -> kCosSolarInclination
  ZEL                      -> kLunarEccentricity                 ZES                      -> kSolarEccentricity
  ZNL                      -> kLunarMeanMotion                   ZNS                      -> kSolarMeanMotion
  ZSINGS                   -> kSinSolarArgPerigee                ZSINI                    -> kSinSolarInclination
  a1                       -> dirCos1                            a10                      -> dirCos10
  a2                       -> dirCos2                            a3                       -> dirCos3
  a4                       -> dirCos4                            a5                       -> dirCos5
  a6                       -> dirCos6                            a7                       -> dirCos7
  a8                       -> dirCos8                            a9                       -> dirCos9
  ainv2                    -> inverseSmaSquared                  aqnv                     -> inverseSemiMajorAxis0
  atime                    -> integratorTime                     betao                    -> sqrtOneMinusEccSq
  betao2                   -> oneMinusEccSq                      bfact                    -> resonanceFrequency
  c                        -> lunarMeanLongitude                 cc                       -> bodyC1
  cosg                     -> cosArgPerigee0                     cosio                    -> cosInclination0
  cosq                     -> cosRaAn0                           ctem                     -> cosNodeLon
  d2201                    -> resonanceD2201                     d2211                    -> resonanceD2211
  d3210                    -> resonanceD3210                     d3222                    -> resonanceD3222
  d4410                    -> resonanceD4410                     d4422                    -> resonanceD4422
  d5220                    -> resonanceD5220                     d5232                    -> resonanceD5232
  d5421                    -> resonanceD5421                     d5433                    -> resonanceD5433
  del1                     -> synchronousDel1                    del2                     -> synchronousDel2
  del3                     -> synchronousDel3                    e3                       -> lunarEcc3
  ee2                      -> lunarEcc2                          eosq                     -> eccentricitySq
  f220                     -> resonanceF220                      f311                     -> resonanceF311
  f330                     -> resonanceF330                      g200                     -> resonanceG200
  g201                     -> resonanceG201                      g211                     -> resonanceG211
  g300                     -> resonanceG300                      g310                     -> resonanceG310
  g322                     -> resonanceG322                      g410                     -> resonanceG410
  g422                     -> resonanceG422                      g520                     -> resonanceG520
  g521                     -> resonanceG521                      g532                     -> resonanceG532
  g533                     -> resonanceG533                      gam                      -> lunarPerigeeLongitude
  gsto                     -> greenwichSiderealTime              jday                     -> daysSince1900
  omgdot                   -> argPerigeeDot                      s1                       -> dot1
  s2                       -> dot2                               s3                       -> dot3
  s4                       -> dot4                               s5                       -> dot5
  s6                       -> dot6                               s7                       -> dot7
  se                       -> bodySecularEcc                     se2                      -> solarEcc2
  se3                      -> solarEcc3                          sgh                      -> bodySecularArgPerigee
  sgh2                     -> solarArgPerigee2                   sgh3                     -> solarArgPerigee3
  sgh4                     -> solarArgPerigee4                   sh2                      -> solarRaAn2
  sh3                      -> solarRaAn3                         shdq                     -> bodySecularRaAn
  si                       -> bodySecularInc                     si2                      -> solarInc2
  si3                      -> solarInc3                          sing                     -> sinArgPerigee0
  sini2                    -> sinSqInc                           sinio                    -> sinInclination0
  sinq                     -> sinRaAn0                           sl                       -> bodySecularLong
  sl2                      -> solarLong2                         sl3                      -> solarLong3
  sl4                      -> solarLong4                         sse                      -> totalSecularEcc
  ssg                      -> totalSecularArgPerigee             ssh                      -> totalSecularRaAn
  ssi                      -> totalSecularInc                    ssl                      -> totalSecularLong
  stem                     -> sinNodeLon                         temp                     -> resonanceCoef
  temp1                    -> resonancePrefactor                 theta2                   -> cosSqInc
  x1                       -> rot1                               x2                       -> rot2
  x3                       -> rot3                               x4                       -> rot4
  x5                       -> rot5                               x6                       -> rot6
  x7                       -> rot7                               x8                       -> rot8
  xfact                    -> resonancePhaseRate                 xgh2                     -> lunarArgPerigee2
  xgh3                     -> lunarArgPerigee3                   xgh4                     -> lunarArgPerigee4
  xh2                      -> lunarRaAn2                         xh3                      -> lunarRaAn3
  xi2                      -> lunarInc2                          xi3                      -> lunarInc3
  xl2                      -> lunarLong2                         xl3                      -> lunarLong3
  xl4                      -> lunarLong4                         xlamo                    -> resonancePhase0
  xli                      -> resonancePhase                     xmdot                    -> meanAnomalyDot
  xni                      -> resonanceMeanMotion                xno2                     -> meanMotionSquared
  xnodce                   -> lunarNodeLongitude                 xnodot                   -> raanDot
  xnoi                     -> inverseMeanMotion                  xpidot                   -> argPerigeeRaAnDot
  z1                       -> grav1                              z11                      -> grav11
  z12                      -> grav12                             z13                      -> grav13
  z2                       -> grav2                              z21                      -> grav21
  z22                      -> grav22                             z23                      -> grav23
  z3                       -> grav3                              z31                      -> grav31
  z32                      -> grav32                             z33                      -> grav33
  zcosg                    -> cosSolarArgPerigee                 zcosgl                   -> cosLunarArgPerigee
  zcosh                    -> cosNodeRef                         zcoshl                   -> cosLunarNodeAngle
  zcosi                    -> cosSolarIncl                       zcosil                   -> cosLunarIncl
  ze                       -> bodyEccentricity                   zmol                     -> lunarMeanAnomaly
  zmos                     -> solarMeanAnomaly                   zn                       -> bodyMeanMotion
  zsing                    -> sinSolarArgPerigee                 zsingl                   -> sinLunarArgPerigee
  zsinh                    -> sinNodeRef                         zsinhl                   -> sinLunarNodeAngle
  zsini                    -> sinSolarIncl                       zsinil                   -> sinLunarIncl
  zx                       -> lunarArgPerigeeSinInterim          zy                       -> lunarArgPerigeeCosInterim
 */
void SGP4::DeepSpaceInitialise(double eccentricitySq,
                               double sinInclination0,
                               double cosInclination0,
                               double sqrtOneMinusEccSq,
                               double cosSqInc,
                               double oneMinusEccSq,
                               double meanAnomalyDot,
                               double argPerigeeDot,
                               double raanDot)
{
    double bodySecularEcc = 0.0;
    double bodySecularInc = 0.0;
    double bodySecularLong = 0.0;
    double bodySecularArgPerigee = 0.0;
    double bodySecularRaAn = 0.0;

    double resonanceFrequency = 0.0;

    const double kSolarMeanMotion = 1.19459E-5;
    const double kSolarC1 = 2.9864797E-6;
    const double kSolarEccentricity = 0.01675;
    const double kLunarMeanMotion = 1.5835218E-4;
    const double kLunarC1 = 4.7968065E-7;
    const double kLunarEccentricity = 0.05490;
    const double kCosSolarInclination = 0.91744867;
    const double kSinSolarInclination = 0.39785416;
    const double kSinSolarArgPerigee = -0.98088458;
    const double kCosSolarArgPerigee = 0.1945905;
    const double kSynchronousQ22 = 1.7891679E-6;
    const double kSynchronousQ31 = 2.1460748E-6;
    const double kSynchronousQ33 = 2.2123015E-7;
    const double kResonanceRoot22 = 1.7891679E-6;
    const double kResonanceRoot32 = 3.7393792E-7;
    const double kResonanceRoot44 = 7.3636953E-9;
    const double kResonanceRoot52 = 1.1428639E-7;
    const double kResonanceRoot54 = 2.1765803E-9;

    const double inverseSemiMajorAxis0 = 1.0 / mElements.RecoveredSemiMajorAxis();
    const double argPerigeeRaAnDot = argPerigeeDot + raanDot;
    const double sinRaAn0 = sin(mElements.AscendingNode());
    const double cosRaAn0 = cos(mElements.AscendingNode());
    const double sinArgPerigee0 = sin(mElements.ArgumentPerigee());
    const double cosArgPerigee0 = cos(mElements.ArgumentPerigee());

    /*
     * initialize lunar / solar terms
     */
    const double daysSince1900 = mElements.Epoch().ToJ1900();

    const double lunarNodeLongitude = Util::WrapTwoPI(4.5236020 - 9.2422029e-4 * daysSince1900);
    const double sinNodeLon = sin(lunarNodeLongitude);
    const double cosNodeLon = cos(lunarNodeLongitude);
    const double cosLunarIncl = 0.91375164 - 0.03568096 * cosNodeLon;
    const double sinLunarIncl = sqrt(1.0 - cosLunarIncl * cosLunarIncl);
    const double sinLunarNodeAngle = 0.089683511 * sinNodeLon / sinLunarIncl;
    const double cosLunarNodeAngle = sqrt(1.0 - sinLunarNodeAngle * sinLunarNodeAngle);
    const double lunarMeanLongitude = 4.7199672 + 0.22997150 * daysSince1900;
    const double lunarPerigeeLongitude = 5.8351514 + 0.0019443680 * daysSince1900;
    mDeepspaceConsts.lunarMeanAnomaly = Util::WrapTwoPI(lunarMeanLongitude - lunarPerigeeLongitude);
    double lunarArgPerigeeSinInterim = 0.39785416 * sinNodeLon / sinLunarIncl;
    double lunarArgPerigeeCosInterim = cosLunarNodeAngle * cosNodeLon + 0.91744867 * sinLunarNodeAngle * sinNodeLon;
    lunarArgPerigeeSinInterim = atan2(lunarArgPerigeeSinInterim, lunarArgPerigeeCosInterim);
    lunarArgPerigeeSinInterim = lunarPerigeeLongitude + lunarArgPerigeeSinInterim - lunarNodeLongitude;

    const double cosLunarArgPerigee = cos(lunarArgPerigeeSinInterim);
    const double sinLunarArgPerigee = sin(lunarArgPerigeeSinInterim);
    mDeepspaceConsts.solarMeanAnomaly = Util::WrapTwoPI(6.2565837 + 0.017201977 * daysSince1900);

    /*
     * do solar terms
     */
    double cosSolarArgPerigee = kCosSolarArgPerigee;
    double sinSolarArgPerigee = kSinSolarArgPerigee;
    double cosSolarIncl = kCosSolarInclination;
    double sinSolarIncl = kSinSolarInclination;
    double cosNodeRef = cosRaAn0;
    double sinNodeRef = sinRaAn0;
    double bodyC1 = kSolarC1;
    double bodyMeanMotion = kSolarMeanMotion;
    double bodyEccentricity = kSolarEccentricity;
    const double inverseMeanMotion = 1.0 / mElements.RecoveredMeanMotion();

    for (int cnt = 0; cnt < 2; cnt++)
    {
        /*
         * solar terms are done a second time after lunar terms are done
         */
        const double dirCos1 = cosSolarArgPerigee * cosNodeRef + sinSolarArgPerigee * cosSolarIncl * sinNodeRef;
        const double dirCos3 = -sinSolarArgPerigee * cosNodeRef + cosSolarArgPerigee * cosSolarIncl * sinNodeRef;
        const double dirCos7 = -cosSolarArgPerigee * sinNodeRef + sinSolarArgPerigee * cosSolarIncl * cosNodeRef;
        const double dirCos8 = sinSolarArgPerigee * sinSolarIncl;
        const double dirCos9 = sinSolarArgPerigee * sinNodeRef + cosSolarArgPerigee * cosSolarIncl * cosNodeRef;
        const double dirCos10 = cosSolarArgPerigee * sinSolarIncl;
        const double dirCos2 = cosInclination0 * dirCos7 + sinInclination0 * dirCos8;
        const double dirCos4 = cosInclination0 * dirCos9 + sinInclination0 * dirCos10;
        const double dirCos5 = -sinInclination0 * dirCos7 + cosInclination0 * dirCos8;
        const double dirCos6 = -sinInclination0 * dirCos9 + cosInclination0 * dirCos10;
        const double rot1 = dirCos1 * cosArgPerigee0 + dirCos2 * sinArgPerigee0;
        const double rot2 = dirCos3 * cosArgPerigee0 + dirCos4 * sinArgPerigee0;
        const double rot3 = -dirCos1 * sinArgPerigee0 + dirCos2 * cosArgPerigee0;
        const double rot4 = -dirCos3 * sinArgPerigee0 + dirCos4 * cosArgPerigee0;
        const double rot5 = dirCos5 * sinArgPerigee0;
        const double rot6 = dirCos6 * sinArgPerigee0;
        const double rot7 = dirCos5 * cosArgPerigee0;
        const double rot8 = dirCos6 * cosArgPerigee0;
        const double grav31 = 12.0 * rot1 * rot1 - 3. * rot3 * rot3;
        const double grav32 = 24.0 * rot1 * rot2 - 6. * rot3 * rot4;
        const double grav33 = 12.0 * rot2 * rot2 - 3. * rot4 * rot4;
        double grav1 = 3.0 * (dirCos1 * dirCos1 + dirCos2 * dirCos2) + grav31 * eccentricitySq;
        double grav2 = 6.0 * (dirCos1 * dirCos3 + dirCos2 * dirCos4) + grav32 * eccentricitySq;
        double grav3 = 3.0 * (dirCos3 * dirCos3 + dirCos4 * dirCos4) + grav33 * eccentricitySq;

        const double grav11 = -6.0 * dirCos1 * dirCos5 + eccentricitySq * (-24. * rot1 * rot7 - 6. * rot3 * rot5);
        const double grav12 = -6.0 * (dirCos1 * dirCos6 + dirCos3 * dirCos5) +
                              eccentricitySq * (-24. * (rot2 * rot7 + rot1 * rot8) - 6. * (rot3 * rot6 + rot4 * rot5));
        const double grav13 = -6.0 * dirCos3 * dirCos6 + eccentricitySq * (-24. * rot2 * rot8 - 6. * rot4 * rot6);
        const double grav21 = 6.0 * dirCos2 * dirCos5 + eccentricitySq * (24. * rot1 * rot5 - 6. * rot3 * rot7);
        const double grav22 = 6.0 * (dirCos4 * dirCos5 + dirCos2 * dirCos6) +
                              eccentricitySq * (24. * (rot2 * rot5 + rot1 * rot6) - 6. * (rot4 * rot7 + rot3 * rot8));
        const double grav23 = 6.0 * dirCos4 * dirCos6 + eccentricitySq * (24. * rot2 * rot6 - 6. * rot4 * rot8);

        grav1 = grav1 + grav1 + oneMinusEccSq * grav31;
        grav2 = grav2 + grav2 + oneMinusEccSq * grav32;
        grav3 = grav3 + grav3 + oneMinusEccSq * grav33;

        const double dot3 = bodyC1 * inverseMeanMotion;
        const double dot2 = -0.5 * dot3 / sqrtOneMinusEccSq;
        const double dot4 = dot3 * sqrtOneMinusEccSq;
        const double dot1 = -15.0 * mElements.Eccentricity() * dot4;
        const double dot5 = rot1 * rot3 + rot2 * rot4;
        const double dot6 = rot2 * rot3 + rot1 * rot4;
        const double dot7 = rot2 * rot4 - rot1 * rot3;

        bodySecularEcc = dot1 * bodyMeanMotion * dot5;
        bodySecularInc = dot2 * bodyMeanMotion * (grav11 + grav13);
        bodySecularLong = -bodyMeanMotion * dot3 * (grav1 + grav3 - 14.0 - 6.0 * eccentricitySq);
        bodySecularArgPerigee = dot4 * bodyMeanMotion * (grav31 + grav33 - 6.0);

        /*
         * replaced
         * sh = -zn * s2 * (z21 + z23
         * with
         * shdq = (-zn * s2 * (z21 + z23)) / sinio
         */
        if (mElements.Inclination() < 5.2359877e-2 || mElements.Inclination() > kPI - 5.2359877e-2)
        {
            bodySecularRaAn = 0.0;
        }
        else
        {
            bodySecularRaAn = (-bodyMeanMotion * dot2 * (grav21 + grav23)) / sinInclination0;
        }

        mDeepspaceConsts.lunarEcc2 = 2.0 * dot1 * dot6;
        mDeepspaceConsts.lunarEcc3 = 2.0 * dot1 * dot7;
        mDeepspaceConsts.lunarInc2 = 2.0 * dot2 * grav12;
        mDeepspaceConsts.lunarInc3 = 2.0 * dot2 * (grav13 - grav11);
        mDeepspaceConsts.lunarLong2 = -2.0 * dot3 * grav2;
        mDeepspaceConsts.lunarLong3 = -2.0 * dot3 * (grav3 - grav1);
        mDeepspaceConsts.lunarLong4 = -2.0 * dot3 * (-21.0 - 9.0 * eccentricitySq) * bodyEccentricity;
        mDeepspaceConsts.lunarArgPerigee2 = 2.0 * dot4 * grav32;
        mDeepspaceConsts.lunarArgPerigee3 = 2.0 * dot4 * (grav33 - grav31);
        mDeepspaceConsts.lunarArgPerigee4 = -18.0 * dot4 * bodyEccentricity;
        mDeepspaceConsts.lunarRaAn2 = -2.0 * dot2 * grav22;
        mDeepspaceConsts.lunarRaAn3 = -2.0 * dot2 * (grav23 - grav21);

        if (cnt == 1)
        {
            break;
        }
        /*
         * do lunar terms
         */
        mDeepspaceConsts.totalSecularEcc = bodySecularEcc;
        mDeepspaceConsts.totalSecularInc = bodySecularInc;
        mDeepspaceConsts.totalSecularLong = bodySecularLong;
        mDeepspaceConsts.totalSecularRaAn = bodySecularRaAn;
        mDeepspaceConsts.totalSecularArgPerigee =
            bodySecularArgPerigee - cosInclination0 * mDeepspaceConsts.totalSecularRaAn;
        mDeepspaceConsts.solarEcc2 = mDeepspaceConsts.lunarEcc2;
        mDeepspaceConsts.solarInc2 = mDeepspaceConsts.lunarInc2;
        mDeepspaceConsts.solarLong2 = mDeepspaceConsts.lunarLong2;
        mDeepspaceConsts.solarArgPerigee2 = mDeepspaceConsts.lunarArgPerigee2;
        mDeepspaceConsts.solarRaAn2 = mDeepspaceConsts.lunarRaAn2;
        mDeepspaceConsts.solarEcc3 = mDeepspaceConsts.lunarEcc3;
        mDeepspaceConsts.solarInc3 = mDeepspaceConsts.lunarInc3;
        mDeepspaceConsts.solarLong3 = mDeepspaceConsts.lunarLong3;
        mDeepspaceConsts.solarArgPerigee3 = mDeepspaceConsts.lunarArgPerigee3;
        mDeepspaceConsts.solarRaAn3 = mDeepspaceConsts.lunarRaAn3;
        mDeepspaceConsts.solarLong4 = mDeepspaceConsts.lunarLong4;
        mDeepspaceConsts.solarArgPerigee4 = mDeepspaceConsts.lunarArgPerigee4;
        cosSolarArgPerigee = cosLunarArgPerigee;
        sinSolarArgPerigee = sinLunarArgPerigee;
        cosSolarIncl = cosLunarIncl;
        sinSolarIncl = sinLunarIncl;
        cosNodeRef = cosLunarNodeAngle * cosRaAn0 + sinLunarNodeAngle * sinRaAn0;
        sinNodeRef = sinRaAn0 * cosLunarNodeAngle - cosRaAn0 * sinLunarNodeAngle;
        bodyMeanMotion = kLunarMeanMotion;
        bodyC1 = kLunarC1;
        bodyEccentricity = kLunarEccentricity;
    }

    mDeepspaceConsts.totalSecularEcc += bodySecularEcc;
    mDeepspaceConsts.totalSecularInc += bodySecularInc;
    mDeepspaceConsts.totalSecularLong += bodySecularLong;
    mDeepspaceConsts.totalSecularArgPerigee += bodySecularArgPerigee - cosInclination0 * bodySecularRaAn;
    mDeepspaceConsts.totalSecularRaAn += bodySecularRaAn;

    mDeepspaceConsts.shape = DeepSpaceConstants::NONE;

    if (mElements.RecoveredMeanMotion() < 0.0052359877 && mElements.RecoveredMeanMotion() > 0.0034906585)
    {
        /*
         * 24h synchronous resonance terms initialisation
         */
        mDeepspaceConsts.shape = DeepSpaceConstants::SYNCHRONOUS;

        const double resonanceG200 = 1.0 + eccentricitySq * (-2.5 + 0.8125 * eccentricitySq);
        const double resonanceG310 = 1.0 + 2.0 * eccentricitySq;
        const double resonanceG300 = 1.0 + eccentricitySq * (-6.0 + 6.60937 * eccentricitySq);
        const double resonanceF220 = 0.75 * (1.0 + cosInclination0) * (1.0 + cosInclination0);
        const double resonanceF311 =
            0.9375 * sinInclination0 * sinInclination0 * (1.0 + 3.0 * cosInclination0) - 0.75 * (1.0 + cosInclination0);
        double resonanceF330 = 1.0 + cosInclination0;
        resonanceF330 = 1.875 * resonanceF330 * resonanceF330 * resonanceF330;
        mDeepspaceConsts.synchronousDel1 = 3.0 * mElements.RecoveredMeanMotion() * mElements.RecoveredMeanMotion() *
                                           inverseSemiMajorAxis0 * inverseSemiMajorAxis0;
        mDeepspaceConsts.synchronousDel2 =
            2.0 * mDeepspaceConsts.synchronousDel1 * resonanceF220 * resonanceG200 * kSynchronousQ22;
        mDeepspaceConsts.synchronousDel3 = 3.0 * mDeepspaceConsts.synchronousDel1 * resonanceF330 * resonanceG300 *
                                           kSynchronousQ33 * inverseSemiMajorAxis0;
        mDeepspaceConsts.synchronousDel1 =
            mDeepspaceConsts.synchronousDel1 * resonanceF311 * resonanceG310 * kSynchronousQ31 * inverseSemiMajorAxis0;

        mDeepspaceConsts.resonancePhase0 =
            Util::WrapTwoPI(mElements.MeanAnomaly() + mElements.AscendingNode() + mElements.ArgumentPerigee() -
                            mDeepspaceConsts.greenwichSiderealTime);
        resonanceFrequency = meanAnomalyDot + argPerigeeRaAnDot - kTHDT + mDeepspaceConsts.totalSecularLong +
                             mDeepspaceConsts.totalSecularArgPerigee + mDeepspaceConsts.totalSecularRaAn;
    }
    else if (mElements.RecoveredMeanMotion() < 8.26e-3 || mElements.RecoveredMeanMotion() > 9.24e-3 ||
             mElements.Eccentricity() < 0.5)
    {
        // do nothing
    }
    else
    {
        /*
         * geopotential resonance initialisation for 12 hour orbits
         */
        mDeepspaceConsts.shape = DeepSpaceConstants::RESONANCE;

        double resonanceG211;
        double resonanceG310;
        double resonanceG322;
        double resonanceG410;
        double resonanceG422;
        double resonanceG520;

        double resonanceG201 = -0.306 - (mElements.Eccentricity() - 0.64) * 0.440;

        if (mElements.Eccentricity() <= 0.65)
        {
            resonanceG211 = EvaluateCubicPolynomial(mElements.Eccentricity(), 3.616, -13.247, 16.290, 0.0);
            resonanceG310 = EvaluateCubicPolynomial(mElements.Eccentricity(), -19.302, 117.390, -228.419, 156.591);
            resonanceG322 = EvaluateCubicPolynomial(mElements.Eccentricity(), -18.9068, 109.7927, -214.6334, 146.5816);
            resonanceG410 = EvaluateCubicPolynomial(mElements.Eccentricity(), -41.122, 242.694, -471.094, 313.953);
            resonanceG422 = EvaluateCubicPolynomial(mElements.Eccentricity(), -146.407, 841.880, -1629.014, 1083.435);
            resonanceG520 = EvaluateCubicPolynomial(mElements.Eccentricity(), -532.114, 3017.977, -5740.032, 3708.276);
        }
        else
        {
            resonanceG211 = EvaluateCubicPolynomial(mElements.Eccentricity(), -72.099, 331.819, -508.738, 266.724);
            resonanceG310 = EvaluateCubicPolynomial(mElements.Eccentricity(), -346.844, 1582.851, -2415.925, 1246.113);
            resonanceG322 = EvaluateCubicPolynomial(mElements.Eccentricity(), -342.585, 1554.908, -2366.899, 1215.972);
            resonanceG410 = EvaluateCubicPolynomial(mElements.Eccentricity(), -1052.797, 4758.686, -7193.992, 3651.957);
            resonanceG422 = EvaluateCubicPolynomial(mElements.Eccentricity(), -3581.69, 16178.11, -24462.77, 12422.52);

            if (mElements.Eccentricity() <= 0.715)
            {
                resonanceG520 = EvaluateCubicPolynomial(mElements.Eccentricity(), 1464.74, -4664.75, 3763.64, 0.0);
            }
            else
            {
                resonanceG520 =
                    EvaluateCubicPolynomial(mElements.Eccentricity(), -5149.66, 29936.92, -54087.36, 31324.56);
            }
        }

        double resonanceG533;
        double resonanceG521;
        double resonanceG532;

        if (mElements.Eccentricity() < 0.7)
        {
            resonanceG533 = EvaluateCubicPolynomial(mElements.Eccentricity(), -919.2277, 4988.61, -9064.77, 5542.21);
            resonanceG521 =
                EvaluateCubicPolynomial(mElements.Eccentricity(), -822.71072, 4568.6173, -8491.4146, 5337.524);
            resonanceG532 = EvaluateCubicPolynomial(mElements.Eccentricity(), -853.666, 4690.25, -8624.77, 5341.4);
        }
        else
        {
            resonanceG533 =
                EvaluateCubicPolynomial(mElements.Eccentricity(), -37995.78, 161616.52, -229838.2, 109377.94);
            resonanceG521 =
                EvaluateCubicPolynomial(mElements.Eccentricity(), -51752.104, 218913.95, -309468.16, 146349.42);
            resonanceG532 =
                EvaluateCubicPolynomial(mElements.Eccentricity(), -40023.88, 170470.89, -242699.48, 115605.82);
        }

        const double sinSqInc = sinInclination0 * sinInclination0;
        const double resonanceF220 = 0.75 * (1.0 + 2.0 * cosInclination0 + cosSqInc);
        const double f221 = 1.5 * sinSqInc;
        const double f321 = 1.875 * sinInclination0 * (1.0 - 2.0 * cosInclination0 - 3.0 * cosSqInc);
        const double f322 = -1.875 * sinInclination0 * (1.0 + 2.0 * cosInclination0 - 3.0 * cosSqInc);
        const double f441 = 35.0 * sinSqInc * resonanceF220;
        const double f442 = 39.3750 * sinSqInc * sinSqInc;
        const double f522 = 9.84375 * sinInclination0 *
                            (sinSqInc * (1.0 - 2.0 * cosInclination0 - 5.0 * cosSqInc) +
                             0.33333333 * (-2.0 + 4.0 * cosInclination0 + 6.0 * cosSqInc));
        const double f523 =
            sinInclination0 * (4.92187512 * sinSqInc * (-2.0 - 4.0 * cosInclination0 + 10.0 * cosSqInc) +
                               6.56250012 * (1.0 + 2.0 * cosInclination0 - 3.0 * cosSqInc));
        const double f542 =
            29.53125 * sinInclination0 *
            (2.0 - 8.0 * cosInclination0 + cosSqInc * (-12.0 + 8.0 * cosInclination0 + 10.0 * cosSqInc));
        const double f543 =
            29.53125 * sinInclination0 *
            (-2.0 - 8.0 * cosInclination0 + cosSqInc * (12.0 + 8.0 * cosInclination0 - 10.0 * cosSqInc));

        const double meanMotionSquared = mElements.RecoveredMeanMotion() * mElements.RecoveredMeanMotion();
        const double inverseSmaSquared = inverseSemiMajorAxis0 * inverseSemiMajorAxis0;

        double resonancePrefactor = 3.0 * meanMotionSquared * inverseSmaSquared;
        double resonanceCoef = resonancePrefactor * kResonanceRoot22;
        mDeepspaceConsts.resonanceD2201 = resonanceCoef * resonanceF220 * resonanceG201;
        mDeepspaceConsts.resonanceD2211 = resonanceCoef * f221 * resonanceG211;

        resonancePrefactor *= inverseSemiMajorAxis0;
        resonanceCoef = resonancePrefactor * kResonanceRoot32;
        mDeepspaceConsts.resonanceD3210 = resonanceCoef * f321 * resonanceG310;
        mDeepspaceConsts.resonanceD3222 = resonanceCoef * f322 * resonanceG322;

        resonancePrefactor *= inverseSemiMajorAxis0;
        resonanceCoef = 2.0 * resonancePrefactor * kResonanceRoot44;
        mDeepspaceConsts.resonanceD4410 = resonanceCoef * f441 * resonanceG410;
        mDeepspaceConsts.resonanceD4422 = resonanceCoef * f442 * resonanceG422;

        resonancePrefactor *= inverseSemiMajorAxis0;
        resonanceCoef = resonancePrefactor * kResonanceRoot52;
        mDeepspaceConsts.resonanceD5220 = resonanceCoef * f522 * resonanceG520;
        mDeepspaceConsts.resonanceD5232 = resonanceCoef * f523 * resonanceG532;

        resonanceCoef = 2.0 * resonancePrefactor * kResonanceRoot54;
        mDeepspaceConsts.resonanceD5421 = resonanceCoef * f542 * resonanceG521;
        mDeepspaceConsts.resonanceD5433 = resonanceCoef * f543 * resonanceG533;

        mDeepspaceConsts.resonancePhase0 =
            Util::WrapTwoPI(mElements.MeanAnomaly() + mElements.AscendingNode() + mElements.AscendingNode() -
                            mDeepspaceConsts.greenwichSiderealTime - mDeepspaceConsts.greenwichSiderealTime);
        resonanceFrequency = meanAnomalyDot + raanDot + raanDot - kTHDT - kTHDT + mDeepspaceConsts.totalSecularLong +
                             mDeepspaceConsts.totalSecularRaAn + mDeepspaceConsts.totalSecularRaAn;
    }

    if (mDeepspaceConsts.shape != DeepSpaceConstants::NONE)
    {
        /*
         * initialise integrator
         */
        mDeepspaceConsts.resonancePhaseRate = resonanceFrequency - mElements.RecoveredMeanMotion();
        mIntegratorParams.integratorTime = 0.0;
        mIntegratorParams.resonanceMeanMotion = mElements.RecoveredMeanMotion();
        mIntegratorParams.resonancePhase = mDeepspaceConsts.resonancePhase0;
    }
}

/**
 * From DeepSpaceConstants, this uses:
 * zmos, se2, se3, si2, si3, sl2, sl3, sl4, sgh2, sgh3, sgh4, sh2, sh3
 * zmol, ee2,  e3, xi2, xi3, xl2, xl3, xl4, xgh2, xgh3, xgh4, xh2, xh3
 */
/*
 * traceability: short identifier -> descriptive name
  ZEL                      -> kLunarEccentricity                 ZES                      -> kSolarEccentricity
  ZNL                      -> kLunarMeanMotion                   ZNS                      -> kSolarMeanMotion
  alfdp                    -> lyddaneAlpha                       betdp                    -> lyddaneBeta
  cosis                    -> cosIncl                            cosok                    -> cosRaAn
  dalf                     -> dAlphaLyddane                      dbet                     -> dBetaLyddane
  dls                      -> longitudePerturb                   dsConstants              -> deepSpaceConstants
  e3                       -> lunarEcc3                          ee2                      -> lunarEcc2
  em                       -> eccentricity                       f2                       -> fourier2
  f3                       -> fourier3                           oldxnodes                -> previousRaAn
  omgasm                   -> argPerigee                         pe                       -> perturbEcc
  pgh                      -> perturbArgPerigee                  ph                       -> perturbRaAn
  pinc                     -> perturbInc                         pl                       -> perturbLong
  se2                      -> solarEcc2                          se3                      -> solarEcc3
  sel                      -> lunarPerturbEcc                    ses                      -> solarPerturbEcc
  sgh2                     -> solarArgPerigee2                   sgh3                     -> solarArgPerigee3
  sgh4                     -> solarArgPerigee4                   sghl                     -> lunarPerturbArgPerigee
  sghs                     -> solarPerturbArgPerigee             sh2                      -> solarRaAn2
  sh3                      -> solarRaAn3                         shl                      -> lunarPerturbRaAn
  shs                      -> solarPerturbRaAn                   si2                      -> solarInc2
  si3                      -> solarInc3                          sil                      -> lunarPerturbInc
  sinis                    -> sinIncl                            sinok                    -> sinRaAn
  sinzf                    -> sinEccentricAnomaly                sis                      -> solarPerturbInc
  sl2                      -> solarLong2                         sl3                      -> solarLong3
  sl4                      -> solarLong4                         sll                      -> lunarPerturbLong
  sls                      -> solarPerturbLong                   xgh2                     -> lunarArgPerigee2
  xgh3                     -> lunarArgPerigee3                   xgh4                     -> lunarArgPerigee4
  xh2                      -> lunarRaAn2                         xh3                      -> lunarRaAn3
  xi2                      -> lunarInc2                          xi3                      -> lunarInc3
  xinc                     -> inclination                        xl2                      -> lunarLong2
  xl3                      -> lunarLong3                         xl4                      -> lunarLong4
  xll                      -> meanAnomaly                        xls                      -> perturbedLongitude
  xnodes                   -> raan                               zf                       -> bodyEccentricAnomaly
  zm                       -> bodyMeanAnomaly                    zmol                     -> lunarMeanAnomaly
  zmos                     -> solarMeanAnomaly
 */
void SGP4::DeepSpacePeriodics(double tsince,
                              const DeepSpaceConstants& deepSpaceConstants,
                              double& eccentricity,
                              double& inclination,
                              double& argPerigee,
                              double& raan,
                              double& meanAnomaly)
{
    const double kSolarEccentricity = 0.01675;
    const double kSolarMeanMotion = 1.19459E-5;
    const double kLunarMeanMotion = 1.5835218E-4;
    const double kLunarEccentricity = 0.05490;

    // calculate solar terms for time tsince
    double bodyMeanAnomaly = Util::WrapTwoPI(deepSpaceConstants.solarMeanAnomaly + kSolarMeanMotion * tsince);
    double bodyEccentricAnomaly = bodyMeanAnomaly + 2.0 * kSolarEccentricity * sin(bodyMeanAnomaly);
    double sinEccentricAnomaly = sin(bodyEccentricAnomaly);
    double fourier2 = 0.5 * sinEccentricAnomaly * sinEccentricAnomaly - 0.25;
    double fourier3 = -0.5 * sinEccentricAnomaly * cos(bodyEccentricAnomaly);

    const double solarPerturbEcc = deepSpaceConstants.solarEcc2 * fourier2 + deepSpaceConstants.solarEcc3 * fourier3;
    const double solarPerturbInc = deepSpaceConstants.solarInc2 * fourier2 + deepSpaceConstants.solarInc3 * fourier3;
    const double solarPerturbLong = deepSpaceConstants.solarLong2 * fourier2 +
                                    deepSpaceConstants.solarLong3 * fourier3 +
                                    deepSpaceConstants.solarLong4 * sinEccentricAnomaly;
    const double solarPerturbArgPerigee = deepSpaceConstants.solarArgPerigee2 * fourier2 +
                                          deepSpaceConstants.solarArgPerigee3 * fourier3 +
                                          deepSpaceConstants.solarArgPerigee4 * sinEccentricAnomaly;
    const double solarPerturbRaAn = deepSpaceConstants.solarRaAn2 * fourier2 + deepSpaceConstants.solarRaAn3 * fourier3;

    // calculate lunar terms for time tsince
    bodyMeanAnomaly = Util::WrapTwoPI(deepSpaceConstants.lunarMeanAnomaly + kLunarMeanMotion * tsince);
    bodyEccentricAnomaly = bodyMeanAnomaly + 2.0 * kLunarEccentricity * sin(bodyMeanAnomaly);
    sinEccentricAnomaly = sin(bodyEccentricAnomaly);
    fourier2 = 0.5 * sinEccentricAnomaly * sinEccentricAnomaly - 0.25;
    fourier3 = -0.5 * sinEccentricAnomaly * cos(bodyEccentricAnomaly);

    const double lunarPerturbEcc = deepSpaceConstants.lunarEcc2 * fourier2 + deepSpaceConstants.lunarEcc3 * fourier3;
    const double lunarPerturbInc = deepSpaceConstants.lunarInc2 * fourier2 + deepSpaceConstants.lunarInc3 * fourier3;
    const double lunarPerturbLong = deepSpaceConstants.lunarLong2 * fourier2 +
                                    deepSpaceConstants.lunarLong3 * fourier3 +
                                    deepSpaceConstants.lunarLong4 * sinEccentricAnomaly;
    const double lunarPerturbArgPerigee = deepSpaceConstants.lunarArgPerigee2 * fourier2 +
                                          deepSpaceConstants.lunarArgPerigee3 * fourier3 +
                                          deepSpaceConstants.lunarArgPerigee4 * sinEccentricAnomaly;
    const double lunarPerturbRaAn = deepSpaceConstants.lunarRaAn2 * fourier2 + deepSpaceConstants.lunarRaAn3 * fourier3;

    // merge calculated values
    const double perturbEcc = solarPerturbEcc + lunarPerturbEcc;
    const double perturbInc = solarPerturbInc + lunarPerturbInc;
    const double perturbLong = solarPerturbLong + lunarPerturbLong;
    const double perturbArgPerigee = solarPerturbArgPerigee + lunarPerturbArgPerigee;
    const double perturbRaAn = solarPerturbRaAn + lunarPerturbRaAn;

    inclination += perturbInc;
    eccentricity += perturbEcc;

    /* Spacetrack report #3 has sin/cos from before perturbations
     * added to xinc (oldxinc), but apparently report # 6 has then
     * from after they are added.
     * use for strn3
     * if (mElements.Inclination() >= 0.2)
     * use for gsfc
     * if (xinc >= 0.2)
     * (moved from start of function)
     */
    const double sinIncl = sin(inclination);
    const double cosIncl = cos(inclination);

    if (inclination >= 0.2)
    {
        // apply periodics directly
        argPerigee += perturbArgPerigee - cosIncl * perturbRaAn / sinIncl;
        raan += perturbRaAn / sinIncl;
        meanAnomaly += perturbLong;
    }
    else
    {
        // apply periodics with lyddane modification
        const double sinRaAn = sin(raan);
        const double cosRaAn = cos(raan);
        double lyddaneAlpha = sinIncl * sinRaAn;
        double lyddaneBeta = sinIncl * cosRaAn;
        const double dAlphaLyddane = perturbRaAn * cosRaAn + perturbInc * cosIncl * sinRaAn;
        const double dBetaLyddane = -perturbRaAn * sinRaAn + perturbInc * cosIncl * cosRaAn;
        lyddaneAlpha += dAlphaLyddane;
        lyddaneBeta += dBetaLyddane;
        raan = Util::WrapTwoPI(raan);
        double perturbedLongitude = meanAnomaly + argPerigee + cosIncl * raan;
        double longitudePerturb = perturbLong + perturbArgPerigee - perturbInc * raan * sinIncl;
        perturbedLongitude += longitudePerturb;
        const double previousRaAn = raan;
        raan = atan2(lyddaneAlpha, lyddaneBeta);
        /**
         * Get perturbed xnodes in to same quadrant as original.
         * RAAN is in the range of 0 to 360 degrees
         * atan2 is in the range of -180 to 180 degrees
         */
        if (std::abs(previousRaAn - raan) > kPI)
        {
            if (raan < previousRaAn)
            {
                raan += kTWOPI;
            }
            else
            {
                raan -= kTWOPI;
            }
        }

        meanAnomaly += perturbLong;
        argPerigee = perturbedLongitude - meanAnomaly - cosIncl * raan;
    }
}

/*
 * traceability: short identifier -> descriptive name
  FASX2                    -> kResonanceFasx2                    FASX4                    -> kResonanceFasx4
  FASX6                    -> kResonanceFasx6                    G22                      -> kResonanceG22
  G32                      -> kResonanceG32                      G44                      -> kResonanceG44
  G52                      -> kResonanceG52                      G54                      -> kResonanceG54
  STEP                     -> kIntegratorStep                    STEP2                    -> kIntegratorStepSquared
  atime                    -> integratorTime                     cConstants               -> commonConstants
  d2201                    -> resonanceD2201                     d2211                    -> resonanceD2211
  d3210                    -> resonanceD3210                     d3222                    -> resonanceD3222
  d4410                    -> resonanceD4410                     d4422                    -> resonanceD4422
  d5220                    -> resonanceD5220                     d5232                    -> resonanceD5232
  d5421                    -> resonanceD5421                     d5433                    -> resonanceD5433
  del1                     -> synchronousDel1                    del2                     -> synchronousDel2
  del3                     -> synchronousDel3                    delt                     -> integrationStep
  dsConstants              -> deepSpaceConstants                 em                       -> eccentricity
  ft                       -> timeRemaining                      gsto                     -> greenwichSiderealTime
  omgasm                   -> argPerigee                         omgdot                   -> argPerigeeDot
  sse                      -> totalSecularEcc                    ssg                      -> totalSecularArgPerigee
  ssh                      -> totalSecularRaAn                   ssi                      -> totalSecularInc
  ssl                      -> totalSecularLong                   theta                    -> greenwichSiderealTime
  x2li                     -> twoResonancePhase                  x2omi                    -> twoArgPerigeeAtTime
  xfact                    -> resonancePhaseRate                 xinc                     -> inclination
  xlTemp                   -> resonancePhaseAtTime               xlamo                    -> resonancePhase0
  xldot                    -> resonancePhaseDot                  xli                      -> resonancePhase
  xll                      -> meanAnomaly                        xn                       -> meanMotion
  xnddt                    -> meanMotionSecondDerivative         xndot                    -> meanMotionDot
  xni                      -> resonanceMeanMotion                xnodes                   -> raan
  xomi                     -> argPerigeeAtTime
 */
void SGP4::DeepSpaceSecular(double tsince,
                            const OrbitalElements& elements,
                            const CommonConstants& commonConstants,
                            const DeepSpaceConstants& deepSpaceConstants,
                            IntegratorParams& integParams,
                            double& meanAnomaly,
                            double& argPerigee,
                            double& raan,
                            double& eccentricity,
                            double& inclination,
                            double& meanMotion)
{
    const double kResonanceG22 = 5.7686396;
    const double kResonanceG32 = 0.95240898;
    const double kResonanceG44 = 1.8014998;
    const double kResonanceG52 = 1.0508330;
    const double kResonanceG54 = 4.4108898;
    const double kResonanceFasx2 = 0.13130908;
    const double kResonanceFasx4 = 2.8843198;
    const double kResonanceFasx6 = 0.37448087;

    const double kIntegratorStep = 720.0;
    const double kIntegratorStepSquared = 259200.0;

    meanAnomaly += deepSpaceConstants.totalSecularLong * tsince;
    argPerigee += deepSpaceConstants.totalSecularArgPerigee * tsince;
    raan += deepSpaceConstants.totalSecularRaAn * tsince;
    eccentricity += deepSpaceConstants.totalSecularEcc * tsince;
    inclination += deepSpaceConstants.totalSecularInc * tsince;

    if (deepSpaceConstants.shape != DeepSpaceConstants::NONE)
    {
        double meanMotionDot = 0.0;
        double meanMotionSecondDerivative = 0.0;
        double resonancePhaseDot = 0.0;
        /*
         * 1st condition (if tsince is less than one time step from epoch)
         * 2nd condition (if atime and
         *     tsince are of opposite signs, so zero crossing required)
         * 3rd condition (if tsince is closer to zero than
         *     atime, only integrate away from zero)
         */
        if (std::abs(tsince) < kIntegratorStep || tsince * integParams.integratorTime <= 0.0 ||
            std::abs(tsince) < std::abs(integParams.integratorTime))
        {
            // restart back at the epoch
            integParams.integratorTime = 0.0;
            // TODO: check
            integParams.resonanceMeanMotion = elements.RecoveredMeanMotion();
            // TODO: check
            integParams.resonancePhase = deepSpaceConstants.resonancePhase0;
        }

        bool running = true;
        while (running)
        {
            // always calculate dot terms ready for integration beginning
            // from the start of the range which is 'atime'
            if (deepSpaceConstants.shape == DeepSpaceConstants::SYNCHRONOUS)
            {
                meanMotionDot =
                    deepSpaceConstants.synchronousDel1 * sin(integParams.resonancePhase - kResonanceFasx2) +
                    deepSpaceConstants.synchronousDel2 * sin(2.0 * (integParams.resonancePhase - kResonanceFasx4)) +
                    deepSpaceConstants.synchronousDel3 * sin(3.0 * (integParams.resonancePhase - kResonanceFasx6));
                meanMotionSecondDerivative =
                    deepSpaceConstants.synchronousDel1 * cos(integParams.resonancePhase - kResonanceFasx2) +
                    2.0 * deepSpaceConstants.synchronousDel2 *
                        cos(2.0 * (integParams.resonancePhase - kResonanceFasx4)) +
                    3.0 * deepSpaceConstants.synchronousDel3 *
                        cos(3.0 * (integParams.resonancePhase - kResonanceFasx6));
            }
            else
            {
                // TODO: check
                const double argPerigeeAtTime =
                    elements.ArgumentPerigee() + commonConstants.argPerigeeDot * integParams.integratorTime;
                const double twoArgPerigeeAtTime = argPerigeeAtTime + argPerigeeAtTime;
                const double twoResonancePhase = integParams.resonancePhase + integParams.resonancePhase;
                meanMotionDot =
                    deepSpaceConstants.resonanceD2201 *
                        sin(twoArgPerigeeAtTime + integParams.resonancePhase - kResonanceG22) +
                    deepSpaceConstants.resonanceD2211 * sin(integParams.resonancePhase - kResonanceG22) +
                    deepSpaceConstants.resonanceD3210 *
                        sin(argPerigeeAtTime + integParams.resonancePhase - kResonanceG32) +
                    deepSpaceConstants.resonanceD3222 *
                        sin(-argPerigeeAtTime + integParams.resonancePhase - kResonanceG32) +
                    deepSpaceConstants.resonanceD4410 * sin(twoArgPerigeeAtTime + twoResonancePhase - kResonanceG44) +
                    deepSpaceConstants.resonanceD4422 * sin(twoResonancePhase - kResonanceG44) +
                    deepSpaceConstants.resonanceD5220 *
                        sin(argPerigeeAtTime + integParams.resonancePhase - kResonanceG52) +
                    deepSpaceConstants.resonanceD5232 *
                        sin(-argPerigeeAtTime + integParams.resonancePhase - kResonanceG52) +
                    deepSpaceConstants.resonanceD5421 * sin(argPerigeeAtTime + twoResonancePhase - kResonanceG54) +
                    deepSpaceConstants.resonanceD5433 * sin(-argPerigeeAtTime + twoResonancePhase - kResonanceG54);
                meanMotionSecondDerivative =
                    deepSpaceConstants.resonanceD2201 *
                        cos(twoArgPerigeeAtTime + integParams.resonancePhase - kResonanceG22) +
                    deepSpaceConstants.resonanceD2211 * cos(integParams.resonancePhase - kResonanceG22) +
                    deepSpaceConstants.resonanceD3210 *
                        cos(argPerigeeAtTime + integParams.resonancePhase - kResonanceG32) +
                    deepSpaceConstants.resonanceD3222 *
                        cos(-argPerigeeAtTime + integParams.resonancePhase - kResonanceG32) +
                    deepSpaceConstants.resonanceD5220 *
                        cos(argPerigeeAtTime + integParams.resonancePhase - kResonanceG52) +
                    deepSpaceConstants.resonanceD5232 *
                        cos(-argPerigeeAtTime + integParams.resonancePhase - kResonanceG52) +
                    2.0 *
                        (deepSpaceConstants.resonanceD4410 *
                             cos(twoArgPerigeeAtTime + twoResonancePhase - kResonanceG44) +
                         deepSpaceConstants.resonanceD4422 * cos(twoResonancePhase - kResonanceG44) +
                         deepSpaceConstants.resonanceD5421 * cos(argPerigeeAtTime + twoResonancePhase - kResonanceG54) +
                         deepSpaceConstants.resonanceD5433 *
                             cos(-argPerigeeAtTime + twoResonancePhase - kResonanceG54));
            }
            resonancePhaseDot = integParams.resonanceMeanMotion + deepSpaceConstants.resonancePhaseRate;
            meanMotionSecondDerivative *= resonancePhaseDot;

            double timeRemaining = tsince - integParams.integratorTime;
            if (std::abs(timeRemaining) >= kIntegratorStep)
            {
                const double integrationStep = (timeRemaining >= 0.0 ? kIntegratorStep : -kIntegratorStep);
                // integrate by a full step ('delt'), updating the cached
                // values for the new 'atime'
                integParams.resonancePhase = integParams.resonancePhase + resonancePhaseDot * integrationStep +
                                             meanMotionDot * kIntegratorStepSquared;
                integParams.resonanceMeanMotion = integParams.resonanceMeanMotion + meanMotionDot * integrationStep +
                                                  meanMotionSecondDerivative * kIntegratorStepSquared;
                integParams.integratorTime += integrationStep;
            }
            else
            {
                // integrate by the difference 'ft' remaining
                meanMotion = integParams.resonanceMeanMotion + meanMotionDot * timeRemaining +
                             meanMotionSecondDerivative * timeRemaining * timeRemaining * 0.5;
                const double resonancePhaseAtTime = integParams.resonancePhase + resonancePhaseDot * timeRemaining +
                                                    meanMotionDot * timeRemaining * timeRemaining * 0.5;

                const double greenwichSiderealTime =
                    Util::WrapTwoPI(deepSpaceConstants.greenwichSiderealTime + tsince * kTHDT);
                if (deepSpaceConstants.shape == DeepSpaceConstants::SYNCHRONOUS)
                {
                    meanAnomaly = resonancePhaseAtTime + greenwichSiderealTime - raan - argPerigee;
                }
                else
                {
                    meanAnomaly = resonancePhaseAtTime + 2.0 * (greenwichSiderealTime - raan);
                }
                running = false;
            }
        }
    }
}

void SGP4::Reset()
{
    mUseSimpleModel = false;
    mUseDeepSpace = false;

    mCommonConsts = {};
    mNearspaceConsts = {};
    mDeepspaceConsts = {};
    mIntegratorParams = {};
}

} // namespace libsgp4
