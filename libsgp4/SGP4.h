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

#include "DecayedException.h"
#include "Eci.h"
#include "OrbitalElements.h"
#include "SatelliteException.h"
#include "Tle.h"

namespace libsgp4
{

/**
 * @mainpage
 *
 * This documents the SGP4 tracking library.
 */

/**
 * @brief The simplified perturbations model 4 propagater.
 */
class SGP4
{
public:
    explicit SGP4(const Tle& tle)
        : mElements(tle)
    {
        Initialise();
    }

    void SetTle(const Tle& tle);
    Eci FindPosition(double tsince) const;
    Eci FindPosition(const DateTime& date) const;

private:
    struct CommonConstants
    {
        double cosInclination0;
        double sinInclination0;
        double eta;
        double t2Coeff;
        double sinSqInc;
        double threeCosSqIncMinus1;
        double sevenCosSqIncMinus1;
        double ayCoeff;
        double xlCoeff;
        double raanDragCoeff;
        double dragCoeff;
        double dragCoeff4;
        double argPerigeeDot;  // secular rate of argPerigee (radians/sec)
        double raanDot;        // secular rate of raan    (radians/sec)
        double meanAnomalyDot; // secular rate of meanAnomaly (radians/sec)
    };

    struct NearSpaceConstants
    {
        double dragCoeff5;
        double argPerigeeDragCoeff;
        double meanAnomalyDragCoeff;
        double deltaMeanAnomaly0;
        double sinMeanAnomaly0;
        double d2Coeff;
        double d3Coeff;
        double d4Coeff;
        double t3Coeff;
        double t4Coeff;
        double t5Coeff;
    };

    struct DeepSpaceConstants
    {
        double greenwichSiderealTime;
        double lunarMeanAnomaly;
        double solarMeanAnomaly;

        /*
         * lunar / solar constants for epoch
         * applied during DeepSpaceSecular()
         */
        double totalSecularEcc;
        double totalSecularInc;
        double totalSecularLong;
        double totalSecularArgPerigee;
        double totalSecularRaAn;
        /*
         * lunar / solar constants
         * used during DeepSpaceCalculateLunarSolarTerms()
         */
        double solarEcc2;
        double solarInc2;
        double solarLong2;
        double solarArgPerigee2;
        double solarRaAn2;
        double solarEcc3;
        double solarInc3;
        double solarLong3;
        double solarArgPerigee3;
        double solarRaAn3;
        double solarLong4;
        double solarArgPerigee4;
        double lunarEcc2;
        double lunarEcc3;
        double lunarInc2;
        double lunarInc3;
        double lunarLong2;
        double lunarLong3;
        double lunarLong4;
        double lunarArgPerigee2;
        double lunarArgPerigee3;
        double lunarArgPerigee4;
        double lunarRaAn2;
        double lunarRaAn3;
        /*
         * used during DeepSpaceCalcDotTerms()
         */
        double resonanceD2201;
        double resonanceD2211;
        double resonanceD3210;
        double resonanceD3222;
        double resonanceD4410;
        double resonanceD4422;
        double resonanceD5220;
        double resonanceD5232;
        double resonanceD5421;
        double resonanceD5433;
        double synchronousDel1;
        double synchronousDel2;
        double synchronousDel3;
        /*
         * integrator constants
         */
        double resonancePhaseRate;
        double resonancePhase0;

        enum TOrbitShape
        {
            NONE,
            RESONANCE,
            SYNCHRONOUS
        } shape;
    };

    struct IntegratorParams
    {
        /*
         * integrator values
         */
        double resonancePhase;
        double resonanceMeanMotion;
        double integratorTime;
    };

    void Initialise();
    static void RecomputeConstants(double inclination,
                                   double& sinInclination0,
                                   double& cosInclination0,
                                   double& threeCosSqIncMinus1,
                                   double& sinSqInc,
                                   double& sevenCosSqIncMinus1,
                                   double& xlCoeff,
                                   double& ayCoeff);
    Eci FindPositionSDP4(double tsince) const;
    Eci FindPositionSGP4(double tsince) const;
    static Eci CalculateFinalPositionVelocity(const DateTime& date,
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
                                              double sinInclination0);
    /**
     * Deep space initialisation
     */
    void DeepSpaceInitialise(double eccentricitySq,
                             double sinInclination0,
                             double cosInclination0,
                             double sqrtOneMinusEccSq,
                             double cosSqInc,
                             double oneMinusEccSq,
                             double meanAnomalyDot,
                             double argPerigeeDot,
                             double raanDot);
    /**
     * Calculate lunar / solar periodics and apply
     */
    static void DeepSpacePeriodics(double tsince,
                                   const DeepSpaceConstants& deepSpaceConstants,
                                   double& eccentricity,
                                   double& inclination,
                                   double& argPerigee,
                                   double& raan,
                                   double& meanAnomaly);
    /**
     * Deep space secular effects
     */
    static void DeepSpaceSecular(double tsince,
                                 const OrbitalElements& elements,
                                 const CommonConstants& commonConstants,
                                 const DeepSpaceConstants& deepSpaceConstants,
                                 IntegratorParams& integParams,
                                 double& meanAnomaly,
                                 double& argPerigee,
                                 double& raan,
                                 double& eccentricity,
                                 double& inclination,
                                 double& meanMotion);

    /**
     * Reset
     */
    void Reset();

    /*
     * the constants used
     */
    struct CommonConstants mCommonConsts;
    struct NearSpaceConstants mNearspaceConsts;
    struct DeepSpaceConstants mDeepspaceConsts;
    mutable struct IntegratorParams mIntegratorParams;

    /*
     * the orbit data
     */
    OrbitalElements mElements;

    /*
     * flags
     */
    bool mUseSimpleModel;
    bool mUseDeepSpace;
};

} // namespace libsgp4
