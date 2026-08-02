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
        double cosio;
        double sinio;
        double eta;
        double t2cof;
        double x1mth2;
        double x3thm1;
        double x7thm1;
        double aycof;
        double xlcof;
        double xnodcf;
        double c1;
        double c4;
        double omgdot; // secular rate of omega (radians/sec)
        double xnodot; // secular rate of xnode (radians/sec)
        double xmdot;  // secular rate of xmo   (radians/sec)
    };

    struct NearSpaceConstants
    {
        double c5;
        double omgcof;
        double xmcof;
        double delmo;
        double sinmo;
        double d2;
        double d3;
        double d4;
        double t3cof;
        double t4cof;
        double t5cof;
    };

    struct DeepSpaceConstants
    {
        double gsto;
        double zmol;
        double zmos;

        /*
         * lunar / solar constants for epoch
         * applied during DeepSpaceSecular()
         */
        double sse;
        double ssi;
        double ssl;
        double ssg;
        double ssh;
        /*
         * lunar / solar constants
         * used during DeepSpaceCalculateLunarSolarTerms()
         */
        double se2;
        double si2;
        double sl2;
        double sgh2;
        double sh2;
        double se3;
        double si3;
        double sl3;
        double sgh3;
        double sh3;
        double sl4;
        double sgh4;
        double ee2;
        double e3;
        double xi2;
        double xi3;
        double xl2;
        double xl3;
        double xl4;
        double xgh2;
        double xgh3;
        double xgh4;
        double xh2;
        double xh3;
        /*
         * used during DeepSpaceCalcDotTerms()
         */
        double d2201;
        double d2211;
        double d3210;
        double d3222;
        double d4410;
        double d4422;
        double d5220;
        double d5232;
        double d5421;
        double d5433;
        double del1;
        double del2;
        double del3;
        /*
         * integrator constants
         */
        double xfact;
        double xlamo;

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
        double xli;
        double xni;
        double atime;
    };

    void Initialise();
    static void RecomputeConstants(double xinc,
                                   double& sinio,
                                   double& cosio,
                                   double& x3thm1,
                                   double& x1mth2,
                                   double& x7thm1,
                                   double& xlcof,
                                   double& aycof);
    Eci FindPositionSDP4(double tsince) const;
    Eci FindPositionSGP4(double tsince) const;
    static Eci CalculateFinalPositionVelocity(
            const DateTime& date,
            double e,
            double a,
            double omega,
            double xl,
            double xnode,
            double xinc,
            double xlcof,
            double aycof,
            double x3thm1,
            double x1mth2,
            double x7thm1,
            double cosio,
            double sinio);
    /**
     * Deep space initialisation
     */
    void DeepSpaceInitialise(
            double eosq,
            double sinio,
            double cosio,
            double betao,
            double theta2,
            double betao2,
            double xmdot,
            double omgdot,
            double xnodot);
    /**
     * Calculate lunar / solar periodics and apply
     */
    static void DeepSpacePeriodics(
            double tsince,
            const DeepSpaceConstants& dsConstants,
            double& em,
            double& xinc,
            double& omgasm,
            double& xnodes,
            double& xll);
    /**
     * Deep space secular effects
     */
    static void DeepSpaceSecular(
            double tsince,
            const OrbitalElements& elements,
            const CommonConstants& cConstants,
            const DeepSpaceConstants& dsConstants,
            IntegratorParams& integParams,
            double& xll,
            double& omgasm,
            double& xnodes,
            double& em,
            double& xinc,
            double& xn);

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
