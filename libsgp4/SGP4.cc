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
                       mCommonConsts.sinio,
                       mCommonConsts.cosio,
                       mCommonConsts.x3thm1,
                       mCommonConsts.x1mth2,
                       mCommonConsts.x7thm1,
                       mCommonConsts.xlcof,
                       mCommonConsts.aycof);

    const double theta2 = mCommonConsts.cosio * mCommonConsts.cosio;
    const double eosq = mElements.Eccentricity() * mElements.Eccentricity();
    const double betao2 = 1.0 - eosq;
    const double betao = sqrt(betao2);

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
    double s4 = kS;
    double qoms24 = kQOMS2T;
    if (mElements.Perigee() < 156.0)
    {
        s4 = mElements.Perigee() - 78.0;
        if (mElements.Perigee() < 98.0)
        {
            s4 = 20.0;
        }
        qoms24 = pow((120.0 - s4) * kAE / kXKMPER, 4.0);
        s4 = s4 / kXKMPER + kAE;
    }

    /*
     * generate constants
     */
    const double pinvsq =
        1.0 / (mElements.RecoveredSemiMajorAxis() * mElements.RecoveredSemiMajorAxis() * betao2 * betao2);
    const double tsi = 1.0 / (mElements.RecoveredSemiMajorAxis() - s4);
    mCommonConsts.eta = mElements.RecoveredSemiMajorAxis() * mElements.Eccentricity() * tsi;
    const double etasq = mCommonConsts.eta * mCommonConsts.eta;
    const double eeta = mElements.Eccentricity() * mCommonConsts.eta;
    const double psisq = std::abs(1.0 - etasq);
    const double coef = qoms24 * pow(tsi, 4.0);
    const double coef1 = coef / pow(psisq, 3.5);
    const double c2 = coef1 * mElements.RecoveredMeanMotion() *
                      (mElements.RecoveredSemiMajorAxis() * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
                       0.75 * kCK2 * tsi / psisq * mCommonConsts.x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    mCommonConsts.c1 = mElements.BStar() * c2;
    mCommonConsts.c4 = 2.0 * mElements.RecoveredMeanMotion() * coef1 * mElements.RecoveredSemiMajorAxis() * betao2 *
                       (mCommonConsts.eta * (2.0 + 0.5 * etasq) + mElements.Eccentricity() * (0.5 + 2.0 * etasq) -
                        2.0 * kCK2 * tsi / (mElements.RecoveredSemiMajorAxis() * psisq) *
                            (-3.0 * mCommonConsts.x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
                             0.75 * mCommonConsts.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) *
                                 cos(2.0 * mElements.ArgumentPerigee())));
    const double theta4 = theta2 * theta2;
    const double temp1 = 3.0 * kCK2 * pinvsq * mElements.RecoveredMeanMotion();
    const double temp2 = temp1 * kCK2 * pinvsq;
    const double temp3 = 1.25 * kCK4 * pinvsq * pinvsq * mElements.RecoveredMeanMotion();
    mCommonConsts.xmdot = mElements.RecoveredMeanMotion() + 0.5 * temp1 * betao * mCommonConsts.x3thm1 +
                          0.0625 * temp2 * betao * (13.0 - 78.0 * theta2 + 137.0 * theta4);
    const double x1m5th = 1.0 - 5.0 * theta2;
    mCommonConsts.omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) +
                           temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
    const double xhdot1 = -temp1 * mCommonConsts.cosio;
    mCommonConsts.xnodot =
        xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 * (3.0 - 7.0 * theta2)) * mCommonConsts.cosio;
    mCommonConsts.xnodcf = 3.5 * betao2 * xhdot1 * mCommonConsts.c1;
    mCommonConsts.t2cof = 1.5 * mCommonConsts.c1;

    if (mUseDeepSpace)
    {
        mDeepspaceConsts.gsto = mElements.Epoch().ToGreenwichSiderealTime();

        DeepSpaceInitialise(eosq,
                            mCommonConsts.sinio,
                            mCommonConsts.cosio,
                            betao,
                            theta2,
                            betao2,
                            mCommonConsts.xmdot,
                            mCommonConsts.omgdot,
                            mCommonConsts.xnodot);
    }
    else
    {
        double c3 = 0.0;
        if (mElements.Eccentricity() > 1.0e-4)
        {
            c3 = coef * tsi * kA3OVK2 * mElements.RecoveredMeanMotion() * kAE * mCommonConsts.sinio /
                 mElements.Eccentricity();
        }

        mNearspaceConsts.c5 =
            2.0 * coef1 * mElements.RecoveredSemiMajorAxis() * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
        mNearspaceConsts.omgcof = mElements.BStar() * c3 * cos(mElements.ArgumentPerigee());

        mNearspaceConsts.xmcof = 0.0;
        if (mElements.Eccentricity() > 1.0e-4)
        {
            mNearspaceConsts.xmcof = -kTWOTHIRD * coef * mElements.BStar() * kAE / eeta;
        }

        mNearspaceConsts.delmo = pow(1.0 + mCommonConsts.eta * (cos(mElements.MeanAnomaly())), 3.0);
        mNearspaceConsts.sinmo = sin(mElements.MeanAnomaly());

        if (!mUseSimpleModel)
        {
            const double c1sq = mCommonConsts.c1 * mCommonConsts.c1;
            mNearspaceConsts.d2 = 4.0 * mElements.RecoveredSemiMajorAxis() * tsi * c1sq;
            const double temp = mNearspaceConsts.d2 * tsi * mCommonConsts.c1 / 3.0;
            mNearspaceConsts.d3 = (17.0 * mElements.RecoveredSemiMajorAxis() + s4) * temp;
            mNearspaceConsts.d4 = 0.5 * temp * mElements.RecoveredSemiMajorAxis() * tsi *
                                  (221.0 * mElements.RecoveredSemiMajorAxis() + 31.0 * s4) * mCommonConsts.c1;
            mNearspaceConsts.t3cof = mNearspaceConsts.d2 + 2.0 * c1sq;
            mNearspaceConsts.t4cof =
                0.25 * (3.0 * mNearspaceConsts.d3 + mCommonConsts.c1 * (12.0 * mNearspaceConsts.d2 + 10.0 * c1sq));
            mNearspaceConsts.t5cof = 0.2 * (3.0 * mNearspaceConsts.d4 + 12.0 * mCommonConsts.c1 * mNearspaceConsts.d3 +
                                            6.0 * mNearspaceConsts.d2 * mNearspaceConsts.d2 +
                                            15.0 * c1sq * (2.0 * mNearspaceConsts.d2 + c1sq));
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

Eci SGP4::FindPositionSDP4(double tsince) const
{
    /*
     * the final values
     */
    double e;
    double a;
    double omega;
    double xl;
    double xnode;
    double xinc;

    /*
     * update for secular gravity and atmospheric drag
     */
    double xmdf = mElements.MeanAnomaly() + mCommonConsts.xmdot * tsince;
    double omgadf = mElements.ArgumentPerigee() + mCommonConsts.omgdot * tsince;
    const double xnoddf = mElements.AscendingNode() + mCommonConsts.xnodot * tsince;

    const double tsq = tsince * tsince;
    xnode = xnoddf + mCommonConsts.xnodcf * tsq;
    double tempa = 1.0 - mCommonConsts.c1 * tsince;
    double tempe = mElements.BStar() * mCommonConsts.c4 * tsince;
    double templ = mCommonConsts.t2cof * tsq;

    double xn = mElements.RecoveredMeanMotion();
    double em = mElements.Eccentricity();
    xinc = mElements.Inclination();

    DeepSpaceSecular(tsince,
                     mElements,
                     mCommonConsts,
                     mDeepspaceConsts,
                     mIntegratorParams,
                     xmdf,
                     omgadf,
                     xnode,
                     em,
                     xinc,
                     xn);

    if (!(xn > 0.0))
    {
        throw SatelliteException("Error: (xn <= 0.0 or not finite)");
    }

    const double v = kXKE / xn;
    a = cbrt(v * v) * tempa * tempa;
    e = em - tempe;
    double xmam = xmdf + mElements.RecoveredMeanMotion() * templ;

    DeepSpacePeriodics(tsince, mDeepspaceConsts, e, xinc, omgadf, xnode, xmam);

    /*
     * keeping xinc positive important unless you need to display xinc
     * and dislike negative inclinations
     */
    if (xinc < 0.0)
    {
        xinc = -xinc;
        xnode += kPI;
        omgadf -= kPI;
    }

    xl = xmam + omgadf + xnode;
    omega = omgadf;

    /*
     * fix tolerance for error recognition
     */
    if (!(e > -0.001))
    {
        throw SatelliteException("Error: (e <= -0.001 or not finite)");
    }
    else if (e < 1.0e-6)
    {
        e = 1.0e-6;
    }
    else if (e > (1.0 - 1.0e-6))
    {
        e = 1.0 - 1.0e-6;
    }

    /*
     * re-compute the perturbed values
     */
    double perturbedSinio;
    double perturbedCosio;
    double perturbedX3thm1;
    double perturbedX1mth2;
    double perturbedX7thm1;
    double perturbedXlcof;
    double perturbedAycof;
    RecomputeConstants(xinc,
                       perturbedSinio,
                       perturbedCosio,
                       perturbedX3thm1,
                       perturbedX1mth2,
                       perturbedX7thm1,
                       perturbedXlcof,
                       perturbedAycof);

    /*
     * using calculated values, find position and velocity
     */
    return CalculateFinalPositionVelocity(mElements.Epoch().AddMinutes(tsince),
                                          e,
                                          a,
                                          omega,
                                          xl,
                                          xnode,
                                          xinc,
                                          perturbedXlcof,
                                          perturbedAycof,
                                          perturbedX3thm1,
                                          perturbedX1mth2,
                                          perturbedX7thm1,
                                          perturbedCosio,
                                          perturbedSinio);
}

void SGP4::RecomputeConstants(double xinc,
                              double& sinio,
                              double& cosio,
                              double& x3thm1,
                              double& x1mth2,
                              double& x7thm1,
                              double& xlcof,
                              double& aycof)
{
    sinio = sin(xinc);
    cosio = cos(xinc);

    const double theta2 = cosio * cosio;

    x3thm1 = 3.0 * theta2 - 1.0;
    x1mth2 = 1.0 - theta2;
    x7thm1 = 7.0 * theta2 - 1.0;

    if (std::abs(cosio + 1.0) > 1.5e-12)
    {
        xlcof = 0.125 * kA3OVK2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    }
    else
    {
        xlcof = 0.125 * kA3OVK2 * sinio * (3.0 + 5.0 * cosio) / 1.5e-12;
    }

    aycof = 0.25 * kA3OVK2 * sinio;
}

Eci SGP4::FindPositionSGP4(double tsince) const
{
    /*
     * the final values
     */
    double e;
    double a;
    double omega;
    double xl;
    double xnode;
    const double xinc = mElements.Inclination();

    /*
     * update for secular gravity and atmospheric drag
     */
    const double xmdf = mElements.MeanAnomaly() + mCommonConsts.xmdot * tsince;
    const double omgadf = mElements.ArgumentPerigee() + mCommonConsts.omgdot * tsince;
    const double xnoddf = mElements.AscendingNode() + mCommonConsts.xnodot * tsince;

    omega = omgadf;
    double xmp = xmdf;

    const double tsq = tsince * tsince;
    xnode = xnoddf + mCommonConsts.xnodcf * tsq;
    double tempa = 1.0 - mCommonConsts.c1 * tsince;
    double tempe = mElements.BStar() * mCommonConsts.c4 * tsince;
    double templ = mCommonConsts.t2cof * tsq;

    if (!mUseSimpleModel)
    {
        const double delomg = mNearspaceConsts.omgcof * tsince;
        const double x1p = 1.0 + mCommonConsts.eta * cos(Util::WrapTwoPI(xmdf));
        const double delm = mNearspaceConsts.xmcof * (x1p * x1p * x1p - mNearspaceConsts.delmo);
        const double temp = delomg + delm;

        xmp += temp;
        omega -= temp;

        const double tcube = tsq * tsince;
        const double tfour = tsince * tcube;

        tempa = tempa - mNearspaceConsts.d2 * tsq - mNearspaceConsts.d3 * tcube - mNearspaceConsts.d4 * tfour;
        tempe += mElements.BStar() * mNearspaceConsts.c5 * (sin(Util::WrapTwoPI(xmp)) - mNearspaceConsts.sinmo);
        templ += mNearspaceConsts.t3cof * tcube + tfour * (mNearspaceConsts.t4cof + tsince * mNearspaceConsts.t5cof);
    }

    a = mElements.RecoveredSemiMajorAxis() * tempa * tempa;
    e = mElements.Eccentricity() - tempe;
    xl = xmp + omega + xnode + mElements.RecoveredMeanMotion() * templ;

    /*
     * fix tolerance for error recognition
     */
    if (!(e > -0.001))
    {
        throw SatelliteException("Error: (e <= -0.001 or not finite)");
    }
    else if (e < 1.0e-6)
    {
        e = 1.0e-6;
    }
    else if (e > (1.0 - 1.0e-6))
    {
        e = 1.0 - 1.0e-6;
    }

    /*
     * using calculated values, find position and velocity
     * we can pass in constants from Initialise() as these dont change
     */
    return CalculateFinalPositionVelocity(mElements.Epoch().AddMinutes(tsince),
                                          e,
                                          a,
                                          omega,
                                          xl,
                                          xnode,
                                          xinc,
                                          mCommonConsts.xlcof,
                                          mCommonConsts.aycof,
                                          mCommonConsts.x3thm1,
                                          mCommonConsts.x1mth2,
                                          mCommonConsts.x7thm1,
                                          mCommonConsts.cosio,
                                          mCommonConsts.sinio);
}

Eci SGP4::CalculateFinalPositionVelocity(const DateTime& dt,
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
                                         double sinio)
{
    if (!std::isfinite(e) || !std::isfinite(a) || !std::isfinite(omega) || !std::isfinite(xl) ||
        !std::isfinite(xnode) || !std::isfinite(xinc))
    {
        throw SatelliteException("Error: (element not finite)");
    }

    const double beta2 = 1.0 - e * e;
    const double sqrtA = sqrt(a);
    const double xn = kXKE / (a * sqrtA);
    /*
     * long period periodics
     */
    const double axn = e * cos(omega);
    const double temp11 = 1.0 / (a * beta2);
    const double xll = temp11 * xlcof * axn;
    const double aynl = temp11 * aycof;
    const double xlt = xl + xll;
    const double ayn = e * sin(omega) + aynl;
    const double elsq = axn * axn + ayn * ayn;

    if (!(elsq < 1.0))
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
    const double capu = fmod(xlt - xnode, kTWOPI);
    double epw = capu;

    double sinepw = 0.0;
    double cosepw = 0.0;
    double ecose = 0.0;
    double esine = 0.0;

    /*
     * sensibility check for N-R correction
     */
    const double maxNewtonNaphson = 1.25 * std::abs(sqrt(elsq));

    bool keplerRunning = true;

    for (int i = 0; i < 10 && keplerRunning; i++)
    {
        sinepw = sin(epw);
        cosepw = cos(epw);
        ecose = axn * cosepw + ayn * sinepw;
        esine = axn * sinepw - ayn * cosepw;

        double f = capu - epw + esine;

        if (std::abs(f) < 1.0e-12)
        {
            keplerRunning = false;
        }
        else
        {
            /*
             * 1st order Newton-Raphson correction
             */
            const double fdot = 1.0 - ecose;
            double deltaEpw = f / fdot;

            /*
             * 2nd order Newton-Raphson correction.
             * f / (fdot - 0.5 * d2f * f/fdot)
             */
            if (i == 0)
            {
                if (deltaEpw > maxNewtonNaphson)
                {
                    deltaEpw = maxNewtonNaphson;
                }
                else if (deltaEpw < -maxNewtonNaphson)
                {
                    deltaEpw = -maxNewtonNaphson;
                }
            }
            else
            {
                deltaEpw = f / (fdot + 0.5 * esine * deltaEpw);
            }

            /*
             * Newton-Raphson correction of -F/DF
             */
            epw += deltaEpw;
        }
    }
    /*
     * short period preliminary quantities
     */
    const double temp21 = 1.0 - elsq;
    const double pl = a * temp21;

    if (!(pl >= 0.0))
    {
        throw SatelliteException("Error: (pl < 0.0 or not finite)");
    }

    const double r = a * (1.0 - ecose);
    const double temp31 = 1.0 / r;
    const double rdot = kXKE * sqrtA * esine * temp31;
    const double rfdot = kXKE * sqrt(pl) * temp31;
    const double temp32 = a * temp31;
    const double betal = sqrt(temp21);
    const double temp33 = 1.0 / (1.0 + betal);
    const double cosu = temp32 * (cosepw - axn + ayn * esine * temp33);
    const double sinu = temp32 * (sinepw - ayn - axn * esine * temp33);
    const double sin2u = 2.0 * sinu * cosu;
    const double cos2u = 2.0 * cosu * cosu - 1.0;

    /*
     * update for short periodics
     */
    const double temp41 = 1.0 / pl;
    const double temp42 = kCK2 * temp41;
    const double temp43 = temp42 * temp41;

    const double rk = r * (1.0 - 1.5 * temp43 * betal * x3thm1) + 0.5 * temp42 * x1mth2 * cos2u;
    const double delta = -0.25 * temp43 * x7thm1 * sin2u;
    const double xnodek = xnode + 1.5 * temp43 * cosio * sin2u;
    const double xinck = xinc + 1.5 * temp43 * cosio * sinio * cos2u;
    const double rdotk = rdot - xn * temp42 * x1mth2 * sin2u;
    const double rfdotk = rfdot + xn * temp42 * (x1mth2 * cos2u + 1.5 * x3thm1);

    /*
     * orientation vectors
     */
    const double sindelta = sin(delta);
    const double cosdelta = cos(delta);
    const double sinuk = sinu * cosdelta + cosu * sindelta;
    const double cosuk = cosu * cosdelta - sinu * sindelta;
    const double sinik = sin(xinck);
    const double cosik = cos(xinck);
    const double sinnok = sin(xnodek);
    const double cosnok = cos(xnodek);
    const double xmx = -sinnok * cosik;
    const double xmy = cosnok * cosik;
    const double ux = xmx * sinuk + cosnok * cosuk;
    const double uy = xmy * sinuk + sinnok * cosuk;
    const double uz = sinik * sinuk;
    const double vx = xmx * cosuk - cosnok * sinuk;
    const double vy = xmy * cosuk - sinnok * sinuk;
    const double vz = sinik * cosuk;
    /*
     * position and velocity
     */
    const double x = rk * ux * kXKMPER;
    const double y = rk * uy * kXKMPER;
    const double z = rk * uz * kXKMPER;
    Vector position(x, y, z);
    const double xdot = (rdotk * ux + rfdotk * vx) * kXKMPER / 60.0;
    const double ydot = (rdotk * uy + rfdotk * vy) * kXKMPER / 60.0;
    const double zdot = (rdotk * uz + rfdotk * vz) * kXKMPER / 60.0;
    Vector velocity(xdot, ydot, zdot);

    if (!(rk >= 1.0))
    {
        if (!std::isfinite(rk))
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

void SGP4::DeepSpaceInitialise(double eosq,
                               double sinio,
                               double cosio,
                               double betao,
                               double theta2,
                               double betao2,
                               double xmdot,
                               double omgdot,
                               double xnodot)
{
    double se = 0.0;
    double si = 0.0;
    double sl = 0.0;
    double sgh = 0.0;
    double shdq = 0.0;

    double bfact = 0.0;

    const double ZNS = 1.19459E-5;
    const double C1SS = 2.9864797E-6;
    const double ZES = 0.01675;
    const double ZNL = 1.5835218E-4;
    const double C1L = 4.7968065E-7;
    const double ZEL = 0.05490;
    const double ZCOSIS = 0.91744867;
    const double ZSINI = 0.39785416;
    const double ZSINGS = -0.98088458;
    const double ZCOSGS = 0.1945905;
    const double Q22 = 1.7891679E-6;
    const double Q31 = 2.1460748E-6;
    const double Q33 = 2.2123015E-7;
    const double ROOT22 = 1.7891679E-6;
    const double ROOT32 = 3.7393792E-7;
    const double ROOT44 = 7.3636953E-9;
    const double ROOT52 = 1.1428639E-7;
    const double ROOT54 = 2.1765803E-9;

    const double aqnv = 1.0 / mElements.RecoveredSemiMajorAxis();
    const double xpidot = omgdot + xnodot;
    const double sinq = sin(mElements.AscendingNode());
    const double cosq = cos(mElements.AscendingNode());
    const double sing = sin(mElements.ArgumentPerigee());
    const double cosg = cos(mElements.ArgumentPerigee());

    /*
     * initialize lunar / solar terms
     */
    const double jday = mElements.Epoch().ToJ1900();

    const double xnodce = Util::WrapTwoPI(4.5236020 - 9.2422029e-4 * jday);
    const double stem = sin(xnodce);
    const double ctem = cos(xnodce);
    const double zcosil = 0.91375164 - 0.03568096 * ctem;
    const double zsinil = sqrt(1.0 - zcosil * zcosil);
    const double zsinhl = 0.089683511 * stem / zsinil;
    const double zcoshl = sqrt(1.0 - zsinhl * zsinhl);
    const double c = 4.7199672 + 0.22997150 * jday;
    const double gam = 5.8351514 + 0.0019443680 * jday;
    mDeepspaceConsts.zmol = Util::WrapTwoPI(c - gam);
    double zx = 0.39785416 * stem / zsinil;
    double zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    zx = atan2(zx, zy);
    zx = gam + zx - xnodce;

    const double zcosgl = cos(zx);
    const double zsingl = sin(zx);
    mDeepspaceConsts.zmos = Util::WrapTwoPI(6.2565837 + 0.017201977 * jday);

    /*
     * do solar terms
     */
    double zcosg = ZCOSGS;
    double zsing = ZSINGS;
    double zcosi = ZCOSIS;
    double zsini = ZSINI;
    double zcosh = cosq;
    double zsinh = sinq;
    double cc = C1SS;
    double zn = ZNS;
    double ze = ZES;
    const double xnoi = 1.0 / mElements.RecoveredMeanMotion();

    for (int cnt = 0; cnt < 2; cnt++)
    {
        /*
         * solar terms are done a second time after lunar terms are done
         */
        const double a1 = zcosg * zcosh + zsing * zcosi * zsinh;
        const double a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
        const double a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
        const double a8 = zsing * zsini;
        const double a9 = zsing * zsinh + zcosg * zcosi * zcosh;
        const double a10 = zcosg * zsini;
        const double a2 = cosio * a7 + sinio * a8;
        const double a4 = cosio * a9 + sinio * a10;
        const double a5 = -sinio * a7 + cosio * a8;
        const double a6 = -sinio * a9 + cosio * a10;
        const double x1 = a1 * cosg + a2 * sing;
        const double x2 = a3 * cosg + a4 * sing;
        const double x3 = -a1 * sing + a2 * cosg;
        const double x4 = -a3 * sing + a4 * cosg;
        const double x5 = a5 * sing;
        const double x6 = a6 * sing;
        const double x7 = a5 * cosg;
        const double x8 = a6 * cosg;
        const double z31 = 12.0 * x1 * x1 - 3. * x3 * x3;
        const double z32 = 24.0 * x1 * x2 - 6. * x3 * x4;
        const double z33 = 12.0 * x2 * x2 - 3. * x4 * x4;
        double z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * eosq;
        double z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * eosq;
        double z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * eosq;

        const double z11 = -6.0 * a1 * a5 + eosq * (-24. * x1 * x7 - 6. * x3 * x5);
        const double z12 = -6.0 * (a1 * a6 + a3 * a5) + eosq * (-24. * (x2 * x7 + x1 * x8) - 6. * (x3 * x6 + x4 * x5));
        const double z13 = -6.0 * a3 * a6 + eosq * (-24. * x2 * x8 - 6. * x4 * x6);
        const double z21 = 6.0 * a2 * a5 + eosq * (24. * x1 * x5 - 6. * x3 * x7);
        const double z22 = 6.0 * (a4 * a5 + a2 * a6) + eosq * (24. * (x2 * x5 + x1 * x6) - 6. * (x4 * x7 + x3 * x8));
        const double z23 = 6.0 * a4 * a6 + eosq * (24. * x2 * x6 - 6. * x4 * x8);

        z1 = z1 + z1 + betao2 * z31;
        z2 = z2 + z2 + betao2 * z32;
        z3 = z3 + z3 + betao2 * z33;

        const double s3 = cc * xnoi;
        const double s2 = -0.5 * s3 / betao;
        const double s4 = s3 * betao;
        const double s1 = -15.0 * mElements.Eccentricity() * s4;
        const double s5 = x1 * x3 + x2 * x4;
        const double s6 = x2 * x3 + x1 * x4;
        const double s7 = x2 * x4 - x1 * x3;

        se = s1 * zn * s5;
        si = s2 * zn * (z11 + z13);
        sl = -zn * s3 * (z1 + z3 - 14.0 - 6.0 * eosq);
        sgh = s4 * zn * (z31 + z33 - 6.0);

        /*
         * replaced
         * sh = -zn * s2 * (z21 + z23
         * with
         * shdq = (-zn * s2 * (z21 + z23)) / sinio
         */
        if (mElements.Inclination() < 5.2359877e-2 || mElements.Inclination() > kPI - 5.2359877e-2)
        {
            shdq = 0.0;
        }
        else
        {
            shdq = (-zn * s2 * (z21 + z23)) / sinio;
        }

        mDeepspaceConsts.ee2 = 2.0 * s1 * s6;
        mDeepspaceConsts.e3 = 2.0 * s1 * s7;
        mDeepspaceConsts.xi2 = 2.0 * s2 * z12;
        mDeepspaceConsts.xi3 = 2.0 * s2 * (z13 - z11);
        mDeepspaceConsts.xl2 = -2.0 * s3 * z2;
        mDeepspaceConsts.xl3 = -2.0 * s3 * (z3 - z1);
        mDeepspaceConsts.xl4 = -2.0 * s3 * (-21.0 - 9.0 * eosq) * ze;
        mDeepspaceConsts.xgh2 = 2.0 * s4 * z32;
        mDeepspaceConsts.xgh3 = 2.0 * s4 * (z33 - z31);
        mDeepspaceConsts.xgh4 = -18.0 * s4 * ze;
        mDeepspaceConsts.xh2 = -2.0 * s2 * z22;
        mDeepspaceConsts.xh3 = -2.0 * s2 * (z23 - z21);

        if (cnt == 1)
        {
            break;
        }
        /*
         * do lunar terms
         */
        mDeepspaceConsts.sse = se;
        mDeepspaceConsts.ssi = si;
        mDeepspaceConsts.ssl = sl;
        mDeepspaceConsts.ssh = shdq;
        mDeepspaceConsts.ssg = sgh - cosio * mDeepspaceConsts.ssh;
        mDeepspaceConsts.se2 = mDeepspaceConsts.ee2;
        mDeepspaceConsts.si2 = mDeepspaceConsts.xi2;
        mDeepspaceConsts.sl2 = mDeepspaceConsts.xl2;
        mDeepspaceConsts.sgh2 = mDeepspaceConsts.xgh2;
        mDeepspaceConsts.sh2 = mDeepspaceConsts.xh2;
        mDeepspaceConsts.se3 = mDeepspaceConsts.e3;
        mDeepspaceConsts.si3 = mDeepspaceConsts.xi3;
        mDeepspaceConsts.sl3 = mDeepspaceConsts.xl3;
        mDeepspaceConsts.sgh3 = mDeepspaceConsts.xgh3;
        mDeepspaceConsts.sh3 = mDeepspaceConsts.xh3;
        mDeepspaceConsts.sl4 = mDeepspaceConsts.xl4;
        mDeepspaceConsts.sgh4 = mDeepspaceConsts.xgh4;
        zcosg = zcosgl;
        zsing = zsingl;
        zcosi = zcosil;
        zsini = zsinil;
        zcosh = zcoshl * cosq + zsinhl * sinq;
        zsinh = sinq * zcoshl - cosq * zsinhl;
        zn = ZNL;
        cc = C1L;
        ze = ZEL;
    }

    mDeepspaceConsts.sse += se;
    mDeepspaceConsts.ssi += si;
    mDeepspaceConsts.ssl += sl;
    mDeepspaceConsts.ssg += sgh - cosio * shdq;
    mDeepspaceConsts.ssh += shdq;

    mDeepspaceConsts.shape = DeepSpaceConstants::NONE;

    if (mElements.RecoveredMeanMotion() < 0.0052359877 && mElements.RecoveredMeanMotion() > 0.0034906585)
    {
        /*
         * 24h synchronous resonance terms initialisation
         */
        mDeepspaceConsts.shape = DeepSpaceConstants::SYNCHRONOUS;

        const double g200 = 1.0 + eosq * (-2.5 + 0.8125 * eosq);
        const double g310 = 1.0 + 2.0 * eosq;
        const double g300 = 1.0 + eosq * (-6.0 + 6.60937 * eosq);
        const double f220 = 0.75 * (1.0 + cosio) * (1.0 + cosio);
        const double f311 = 0.9375 * sinio * sinio * (1.0 + 3.0 * cosio) - 0.75 * (1.0 + cosio);
        double f330 = 1.0 + cosio;
        f330 = 1.875 * f330 * f330 * f330;
        mDeepspaceConsts.del1 = 3.0 * mElements.RecoveredMeanMotion() * mElements.RecoveredMeanMotion() * aqnv * aqnv;
        mDeepspaceConsts.del2 = 2.0 * mDeepspaceConsts.del1 * f220 * g200 * Q22;
        mDeepspaceConsts.del3 = 3.0 * mDeepspaceConsts.del1 * f330 * g300 * Q33 * aqnv;
        mDeepspaceConsts.del1 = mDeepspaceConsts.del1 * f311 * g310 * Q31 * aqnv;

        mDeepspaceConsts.xlamo = Util::WrapTwoPI(mElements.MeanAnomaly() + mElements.AscendingNode() +
                                                 mElements.ArgumentPerigee() - mDeepspaceConsts.gsto);
        bfact = xmdot + xpidot - kTHDT + mDeepspaceConsts.ssl + mDeepspaceConsts.ssg + mDeepspaceConsts.ssh;
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

        double g211;
        double g310;
        double g322;
        double g410;
        double g422;
        double g520;

        double g201 = -0.306 - (mElements.Eccentricity() - 0.64) * 0.440;

        if (mElements.Eccentricity() <= 0.65)
        {
            g211 = EvaluateCubicPolynomial(mElements.Eccentricity(), 3.616, -13.247, 16.290, 0.0);
            g310 = EvaluateCubicPolynomial(mElements.Eccentricity(), -19.302, 117.390, -228.419, 156.591);
            g322 = EvaluateCubicPolynomial(mElements.Eccentricity(), -18.9068, 109.7927, -214.6334, 146.5816);
            g410 = EvaluateCubicPolynomial(mElements.Eccentricity(), -41.122, 242.694, -471.094, 313.953);
            g422 = EvaluateCubicPolynomial(mElements.Eccentricity(), -146.407, 841.880, -1629.014, 1083.435);
            g520 = EvaluateCubicPolynomial(mElements.Eccentricity(), -532.114, 3017.977, -5740.032, 3708.276);
        }
        else
        {
            g211 = EvaluateCubicPolynomial(mElements.Eccentricity(), -72.099, 331.819, -508.738, 266.724);
            g310 = EvaluateCubicPolynomial(mElements.Eccentricity(), -346.844, 1582.851, -2415.925, 1246.113);
            g322 = EvaluateCubicPolynomial(mElements.Eccentricity(), -342.585, 1554.908, -2366.899, 1215.972);
            g410 = EvaluateCubicPolynomial(mElements.Eccentricity(), -1052.797, 4758.686, -7193.992, 3651.957);
            g422 = EvaluateCubicPolynomial(mElements.Eccentricity(), -3581.69, 16178.11, -24462.77, 12422.52);

            if (mElements.Eccentricity() <= 0.715)
            {
                g520 = EvaluateCubicPolynomial(mElements.Eccentricity(), 1464.74, -4664.75, 3763.64, 0.0);
            }
            else
            {
                g520 = EvaluateCubicPolynomial(mElements.Eccentricity(), -5149.66, 29936.92, -54087.36, 31324.56);
            }
        }

        double g533;
        double g521;
        double g532;

        if (mElements.Eccentricity() < 0.7)
        {
            g533 = EvaluateCubicPolynomial(mElements.Eccentricity(), -919.2277, 4988.61, -9064.77, 5542.21);
            g521 = EvaluateCubicPolynomial(mElements.Eccentricity(), -822.71072, 4568.6173, -8491.4146, 5337.524);
            g532 = EvaluateCubicPolynomial(mElements.Eccentricity(), -853.666, 4690.25, -8624.77, 5341.4);
        }
        else
        {
            g533 = EvaluateCubicPolynomial(mElements.Eccentricity(), -37995.78, 161616.52, -229838.2, 109377.94);
            g521 = EvaluateCubicPolynomial(mElements.Eccentricity(), -51752.104, 218913.95, -309468.16, 146349.42);
            g532 = EvaluateCubicPolynomial(mElements.Eccentricity(), -40023.88, 170470.89, -242699.48, 115605.82);
        }

        const double sini2 = sinio * sinio;
        const double f220 = 0.75 * (1.0 + 2.0 * cosio + theta2);
        const double f221 = 1.5 * sini2;
        const double f321 = 1.875 * sinio * (1.0 - 2.0 * cosio - 3.0 * theta2);
        const double f322 = -1.875 * sinio * (1.0 + 2.0 * cosio - 3.0 * theta2);
        const double f441 = 35.0 * sini2 * f220;
        const double f442 = 39.3750 * sini2 * sini2;
        const double f522 =
            9.84375 * sinio *
            (sini2 * (1.0 - 2.0 * cosio - 5.0 * theta2) + 0.33333333 * (-2.0 + 4.0 * cosio + 6.0 * theta2));
        const double f523 = sinio * (4.92187512 * sini2 * (-2.0 - 4.0 * cosio + 10.0 * theta2) +
                                     6.56250012 * (1.0 + 2.0 * cosio - 3.0 * theta2));
        const double f542 = 29.53125 * sinio * (2.0 - 8.0 * cosio + theta2 * (-12.0 + 8.0 * cosio + 10.0 * theta2));
        const double f543 = 29.53125 * sinio * (-2.0 - 8.0 * cosio + theta2 * (12.0 + 8.0 * cosio - 10.0 * theta2));

        const double xno2 = mElements.RecoveredMeanMotion() * mElements.RecoveredMeanMotion();
        const double ainv2 = aqnv * aqnv;

        double temp1 = 3.0 * xno2 * ainv2;
        double temp = temp1 * ROOT22;
        mDeepspaceConsts.d2201 = temp * f220 * g201;
        mDeepspaceConsts.d2211 = temp * f221 * g211;

        temp1 *= aqnv;
        temp = temp1 * ROOT32;
        mDeepspaceConsts.d3210 = temp * f321 * g310;
        mDeepspaceConsts.d3222 = temp * f322 * g322;

        temp1 *= aqnv;
        temp = 2.0 * temp1 * ROOT44;
        mDeepspaceConsts.d4410 = temp * f441 * g410;
        mDeepspaceConsts.d4422 = temp * f442 * g422;

        temp1 *= aqnv;
        temp = temp1 * ROOT52;
        mDeepspaceConsts.d5220 = temp * f522 * g520;
        mDeepspaceConsts.d5232 = temp * f523 * g532;

        temp = 2.0 * temp1 * ROOT54;
        mDeepspaceConsts.d5421 = temp * f542 * g521;
        mDeepspaceConsts.d5433 = temp * f543 * g533;

        mDeepspaceConsts.xlamo =
            Util::WrapTwoPI(mElements.MeanAnomaly() + mElements.AscendingNode() + mElements.AscendingNode() -
                            mDeepspaceConsts.gsto - mDeepspaceConsts.gsto);
        bfact = xmdot + xnodot + xnodot - kTHDT - kTHDT + mDeepspaceConsts.ssl + mDeepspaceConsts.ssh +
                mDeepspaceConsts.ssh;
    }

    if (mDeepspaceConsts.shape != DeepSpaceConstants::NONE)
    {
        /*
         * initialise integrator
         */
        mDeepspaceConsts.xfact = bfact - mElements.RecoveredMeanMotion();
        mIntegratorParams.atime = 0.0;
        mIntegratorParams.xni = mElements.RecoveredMeanMotion();
        mIntegratorParams.xli = mDeepspaceConsts.xlamo;
    }
}

/**
 * From DeepSpaceConstants, this uses:
 * zmos, se2, se3, si2, si3, sl2, sl3, sl4, sgh2, sgh3, sgh4, sh2, sh3
 * zmol, ee2,  e3, xi2, xi3, xl2, xl3, xl4, xgh2, xgh3, xgh4, xh2, xh3
 */
void SGP4::DeepSpacePeriodics(double tsince,
                              const DeepSpaceConstants& dsConstants,
                              double& em,
                              double& xinc,
                              double& omgasm,
                              double& xnodes,
                              double& xll)
{
    const double ZES = 0.01675;
    const double ZNS = 1.19459E-5;
    const double ZNL = 1.5835218E-4;
    const double ZEL = 0.05490;

    // calculate solar terms for time tsince
    double zm = Util::WrapTwoPI(dsConstants.zmos + ZNS * tsince);
    double zf = zm + 2.0 * ZES * sin(zm);
    double sinzf = sin(zf);
    double f2 = 0.5 * sinzf * sinzf - 0.25;
    double f3 = -0.5 * sinzf * cos(zf);

    const double ses = dsConstants.se2 * f2 + dsConstants.se3 * f3;
    const double sis = dsConstants.si2 * f2 + dsConstants.si3 * f3;
    const double sls = dsConstants.sl2 * f2 + dsConstants.sl3 * f3 + dsConstants.sl4 * sinzf;
    const double sghs = dsConstants.sgh2 * f2 + dsConstants.sgh3 * f3 + dsConstants.sgh4 * sinzf;
    const double shs = dsConstants.sh2 * f2 + dsConstants.sh3 * f3;

    // calculate lunar terms for time tsince
    zm = Util::WrapTwoPI(dsConstants.zmol + ZNL * tsince);
    zf = zm + 2.0 * ZEL * sin(zm);
    sinzf = sin(zf);
    f2 = 0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * cos(zf);

    const double sel = dsConstants.ee2 * f2 + dsConstants.e3 * f3;
    const double sil = dsConstants.xi2 * f2 + dsConstants.xi3 * f3;
    const double sll = dsConstants.xl2 * f2 + dsConstants.xl3 * f3 + dsConstants.xl4 * sinzf;
    const double sghl = dsConstants.xgh2 * f2 + dsConstants.xgh3 * f3 + dsConstants.xgh4 * sinzf;
    const double shl = dsConstants.xh2 * f2 + dsConstants.xh3 * f3;

    // merge calculated values
    const double pe = ses + sel;
    const double pinc = sis + sil;
    const double pl = sls + sll;
    const double pgh = sghs + sghl;
    const double ph = shs + shl;

    xinc += pinc;
    em += pe;

    /* Spacetrack report #3 has sin/cos from before perturbations
     * added to xinc (oldxinc), but apparently report # 6 has then
     * from after they are added.
     * use for strn3
     * if (mElements.Inclination() >= 0.2)
     * use for gsfc
     * if (xinc >= 0.2)
     * (moved from start of function)
     */
    const double sinis = sin(xinc);
    const double cosis = cos(xinc);

    if (xinc >= 0.2)
    {
        // apply periodics directly
        omgasm += pgh - cosis * ph / sinis;
        xnodes += ph / sinis;
        xll += pl;
    }
    else
    {
        // apply periodics with lyddane modification
        const double sinok = sin(xnodes);
        const double cosok = cos(xnodes);
        double alfdp = sinis * sinok;
        double betdp = sinis * cosok;
        const double dalf = ph * cosok + pinc * cosis * sinok;
        const double dbet = -ph * sinok + pinc * cosis * cosok;
        alfdp += dalf;
        betdp += dbet;
        xnodes = Util::WrapTwoPI(xnodes);
        double xls = xll + omgasm + cosis * xnodes;
        double dls = pl + pgh - pinc * xnodes * sinis;
        xls += dls;
        const double oldxnodes = xnodes;
        xnodes = atan2(alfdp, betdp);
        /**
         * Get perturbed xnodes in to same quadrant as original.
         * RAAN is in the range of 0 to 360 degrees
         * atan2 is in the range of -180 to 180 degrees
         */
        if (std::abs(oldxnodes - xnodes) > kPI)
        {
            if (xnodes < oldxnodes)
            {
                xnodes += kTWOPI;
            }
            else
            {
                xnodes -= kTWOPI;
            }
        }

        xll += pl;
        omgasm = xls - xll - cosis * xnodes;
    }
}

void SGP4::DeepSpaceSecular(double tsince,
                            const OrbitalElements& elements,
                            const CommonConstants& cConstants,
                            const DeepSpaceConstants& dsConstants,
                            IntegratorParams& integParams,
                            double& xll,
                            double& omgasm,
                            double& xnodes,
                            double& em,
                            double& xinc,
                            double& xn)
{
    const double G22 = 5.7686396;
    const double G32 = 0.95240898;
    const double G44 = 1.8014998;
    const double G52 = 1.0508330;
    const double G54 = 4.4108898;
    const double FASX2 = 0.13130908;
    const double FASX4 = 2.8843198;
    const double FASX6 = 0.37448087;

    const double STEP = 720.0;
    const double STEP2 = 259200.0;

    xll += dsConstants.ssl * tsince;
    omgasm += dsConstants.ssg * tsince;
    xnodes += dsConstants.ssh * tsince;
    em += dsConstants.sse * tsince;
    xinc += dsConstants.ssi * tsince;

    if (dsConstants.shape != DeepSpaceConstants::NONE)
    {
        double xndot = 0.0;
        double xnddt = 0.0;
        double xldot = 0.0;
        /*
         * 1st condition (if tsince is less than one time step from epoch)
         * 2nd condition (if atime and
         *     tsince are of opposite signs, so zero crossing required)
         * 3rd condition (if tsince is closer to zero than
         *     atime, only integrate away from zero)
         */
        if (std::abs(tsince) < STEP || tsince * integParams.atime <= 0.0 ||
            std::abs(tsince) < std::abs(integParams.atime))
        {
            // restart back at the epoch
            integParams.atime = 0.0;
            // TODO: check
            integParams.xni = elements.RecoveredMeanMotion();
            // TODO: check
            integParams.xli = dsConstants.xlamo;
        }

        bool running = true;
        while (running)
        {
            // always calculate dot terms ready for integration beginning
            // from the start of the range which is 'atime'
            if (dsConstants.shape == DeepSpaceConstants::SYNCHRONOUS)
            {
                xndot = dsConstants.del1 * sin(integParams.xli - FASX2) +
                        dsConstants.del2 * sin(2.0 * (integParams.xli - FASX4)) +
                        dsConstants.del3 * sin(3.0 * (integParams.xli - FASX6));
                xnddt = dsConstants.del1 * cos(integParams.xli - FASX2) +
                        2.0 * dsConstants.del2 * cos(2.0 * (integParams.xli - FASX4)) +
                        3.0 * dsConstants.del3 * cos(3.0 * (integParams.xli - FASX6));
            }
            else
            {
                // TODO: check
                const double xomi = elements.ArgumentPerigee() + cConstants.omgdot * integParams.atime;
                const double x2omi = xomi + xomi;
                const double x2li = integParams.xli + integParams.xli;
                xndot = dsConstants.d2201 * sin(x2omi + integParams.xli - G22) +
                        dsConstants.d2211 * sin(integParams.xli - G22) +
                        dsConstants.d3210 * sin(xomi + integParams.xli - G32) +
                        dsConstants.d3222 * sin(-xomi + integParams.xli - G32) +
                        dsConstants.d4410 * sin(x2omi + x2li - G44) + dsConstants.d4422 * sin(x2li - G44) +
                        dsConstants.d5220 * sin(xomi + integParams.xli - G52) +
                        dsConstants.d5232 * sin(-xomi + integParams.xli - G52) +
                        dsConstants.d5421 * sin(xomi + x2li - G54) + dsConstants.d5433 * sin(-xomi + x2li - G54);
                xnddt =
                    dsConstants.d2201 * cos(x2omi + integParams.xli - G22) +
                    dsConstants.d2211 * cos(integParams.xli - G22) +
                    dsConstants.d3210 * cos(xomi + integParams.xli - G32) +
                    dsConstants.d3222 * cos(-xomi + integParams.xli - G32) +
                    dsConstants.d5220 * cos(xomi + integParams.xli - G52) +
                    dsConstants.d5232 * cos(-xomi + integParams.xli - G52) +
                    2.0 * (dsConstants.d4410 * cos(x2omi + x2li - G44) + dsConstants.d4422 * cos(x2li - G44) +
                           dsConstants.d5421 * cos(xomi + x2li - G54) + dsConstants.d5433 * cos(-xomi + x2li - G54));
            }
            xldot = integParams.xni + dsConstants.xfact;
            xnddt *= xldot;

            double ft = tsince - integParams.atime;
            if (std::abs(ft) >= STEP)
            {
                const double delt = (ft >= 0.0 ? STEP : -STEP);
                // integrate by a full step ('delt'), updating the cached
                // values for the new 'atime'
                integParams.xli = integParams.xli + xldot * delt + xndot * STEP2;
                integParams.xni = integParams.xni + xndot * delt + xnddt * STEP2;
                integParams.atime += delt;
            }
            else
            {
                // integrate by the difference 'ft' remaining
                xn = integParams.xni + xndot * ft + xnddt * ft * ft * 0.5;
                const double xlTemp = integParams.xli + xldot * ft + xndot * ft * ft * 0.5;

                const double theta = Util::WrapTwoPI(dsConstants.gsto + tsince * kTHDT);
                if (dsConstants.shape == DeepSpaceConstants::SYNCHRONOUS)
                {
                    xll = xlTemp + theta - xnodes - omgasm;
                }
                else
                {
                    xll = xlTemp + 2.0 * (theta - xnodes);
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
