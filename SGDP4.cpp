#include "SGDP4.h"

#include "Vector.h"
#include "SatelliteException.h"

#include <cmath>
#include <iomanip>

SGDP4::SGDP4(void) {
    first_run_ = true;
    SetConstants(CONSTANTS_WGS72);
}

SGDP4::~SGDP4(void) {
}

void SGDP4::SetConstants(EnumConstants constants) {
    constants_.AE = 1.0;
    constants_.TWOTHRD = 2.0 / 3.0;
    switch (constants) {
        case CONSTANTS_WGS72_OLD:
            constants_.MU = 0.0;
            constants_.XKMPER = 6378.135;
            constants_.XKE = 0.0743669161;
            constants_.XJ2 = 0.001082616;
            constants_.XJ3 = -0.00000253881;
            constants_.XJ4 = -0.00000165597;
            constants_.J3OJ2 = constants_.XJ3 / constants_.XJ2;
            break;
        case CONSTANTS_WGS72:
            constants_.MU = 398600.8;
            constants_.XKMPER = 6378.135;
            constants_.XKE = 60.0 / sqrt(constants_.XKMPER * constants_.XKMPER * constants_.XKMPER / constants_.MU);
            constants_.XJ2 = 0.001082616;
            constants_.XJ3 = -0.00000253881;
            constants_.XJ4 = -0.00000165597;
            constants_.J3OJ2 = constants_.XJ3 / constants_.XJ2;
            break;
        case CONSTANTS_WGS84:
            constants_.MU = 398600.5;
            constants_.XKMPER = 6378.137;
            constants_.XKE = 60.0 / sqrt(constants_.XKMPER * constants_.XKMPER * constants_.XKMPER / constants_.MU);
            constants_.XJ2 = 0.00108262998905;
            constants_.XJ3 = -0.00000253215306;
            constants_.XJ4 = -0.00000161098761;
            constants_.J3OJ2 = constants_.XJ3 / constants_.XJ2;
            break;
        default:
            throw new SatelliteException("Unrecognised constant value");
            break;
    }
    constants_.S = constants_.AE + 78.0 / constants_.XKMPER;
    constants_.CK2 = 0.5 * constants_.XJ2 * constants_.AE * constants_.AE;
    constants_.CK4 = -0.375 * constants_.XJ4 * constants_.AE * constants_.AE * constants_.AE * constants_.AE;
    constants_.QOMS2T = pow((120.0 - 78.0) * constants_.AE / constants_.XKMPER, 4.0);
}

void SGDP4::SetTle(const Tle& tle) {

    /*
     * extract and format tle data
     */
    mean_anomoly_ = tle.GetField(Tle::FLD_M, Tle::U_RAD);
    ascending_node_ = tle.GetField(Tle::FLD_RAAN, Tle::U_RAD);
    argument_perigee_ = tle.GetField(Tle::FLD_ARGPER, Tle::U_RAD);
    eccentricity_ = tle.GetField(Tle::FLD_E);
    inclination_ = tle.GetField(Tle::FLD_I, Tle::U_RAD);
    mean_motion_ = tle.GetField(Tle::FLD_MMOTION) * Globals::TWOPI() / Globals::MIN_PER_DAY();
    bstar_ = tle.GetField(Tle::FLD_BSTAR);
    epoch_ = tle.GetEpoch();

    /*
     * error checks
     */
    if (eccentricity_ < 0.0 || eccentricity_ > 1.0 - 1.0e-3) {
        throw new SatelliteException("Eccentricity out of range");
    }

    if (inclination_ < 0.0 || eccentricity_ > Globals::PI()) {
        throw new SatelliteException("Inclination out of range");
    }

    /*
     * recover original mean motion (xnodp) and semimajor axis (aodp)
     * from input elements
     */
    double a1 = pow(constants_.XKE / MeanMotion(), constants_.TWOTHRD);
    i_cosio_ = cos(Inclination());
    i_sinio_ = sin(Inclination());
    double theta2 = i_cosio_ * i_cosio_;
    i_x3thm1_ = 3.0 * theta2 - 1.0;
    double eosq = Eccentricity() * Eccentricity();
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);
    double temp = (1.5 * constants_.CK2) * i_x3thm1_ / (betao * betao2);
    double del1 = temp / (a1 * a1);
    double a0 = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + del1 * 134.0 / 81.0)));
    double del0 = temp / (a0 * a0);

    recovered_mean_motion_ = MeanMotion() / (1.0 + del0);
    recovered_semi_major_axis_ = a0 / (1.0 - del0);

    /*
     * find perigee and period
     */
    perigee_ = (RecoveredSemiMajorAxis() * (1.0 - Eccentricity()) - constants_.AE) * constants_.XKMPER;
    period_ = Globals::TWOPI() / RecoveredMeanMotion();

    Initialize(theta2, betao2, betao, eosq);
}

void SGDP4::Initialize(const double& theta2, const double& betao2, const double& betao, const double& eosq) {

    if (Period() >= 225.0) {
        i_use_deep_space_ = true;
    } else {
        i_use_deep_space_ = false;
        i_use_simple_model_ = false;
        /*
         * for perigee less than 220 kilometers, the simple_model flag is set and
         * the equations are truncated to linear variation in sqrt a and
         * quadratic variation in mean anomly. also, the c3 term, the
         * delta omega term and the delta m term are dropped
         */
        if (Perigee() < 220.0) {
            i_use_simple_model_ = true;
        }
    }

    /*
     * for perigee below 156km, the values of
     * s4 and qoms2t are altered
     */
    double s4 = constants_.S;
    double qoms24 = constants_.QOMS2T;
    if (Perigee() < 156.0) {
        s4 = Perigee() - 78.0;
        if (Perigee() < 98.0) {
            s4 = 20.0;
        }
        qoms24 = pow((120.0 - s4) * constants_.AE / constants_.XKMPER, 4.0);
        s4 = s4 / constants_.XKMPER + constants_.AE;
    }

    /*
     * generate constants
     */
    double pinvsq = 1.0 / (RecoveredSemiMajorAxis() * RecoveredSemiMajorAxis() * betao2 * betao2);
    double tsi = 1.0 / (RecoveredSemiMajorAxis() - s4);
    i_eta_ = RecoveredSemiMajorAxis() * Eccentricity() * tsi;
    double etasq = i_eta_ * i_eta_;
    double eeta = Eccentricity() * i_eta_;
    double psisq = fabs(1.0 - etasq);
    double coef = qoms24 * pow(tsi, 4.0);
    double coef1 = coef / pow(psisq, 3.5);
    double c2 = coef1 * RecoveredMeanMotion() * (RecoveredSemiMajorAxis() *
            (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
            0.75 * constants_.CK2 * tsi / psisq *
            i_x3thm1_ * (8.0 + 3.0 * etasq *
            (8.0 + etasq)));
    i_c1_ = BStar() * c2;
    i_a3ovk2_ = -constants_.XJ3 / constants_.CK2 * pow(constants_.AE, 3.0);
    i_x1mth2_ = 1.0 - theta2;
    i_c4_ = 2.0 * RecoveredMeanMotion() * coef1 * RecoveredSemiMajorAxis() * betao2 *
            (i_eta_ * (2.0 + 0.5 * etasq) + Eccentricity() * (0.5 + 2.0 * etasq) -
            2.0 * constants_.CK2 * tsi / (RecoveredSemiMajorAxis() * psisq) *
            (-3.0 * i_x3thm1_ * (1.0 - 2.0 * eeta + etasq *
            (1.5 - 0.5 * eeta)) + 0.75 * i_x1mth2_ * (2.0 * etasq - eeta *
            (1.0 + etasq)) * cos(2.0 * ArgumentPerigee())));
    double theta4 = theta2 * theta2;
    double temp1 = 3.0 * constants_.CK2 * pinvsq * RecoveredMeanMotion();
    double temp2 = temp1 * constants_.CK2 * pinvsq;
    double temp3 = 1.25 * constants_.CK4 * pinvsq * pinvsq * RecoveredMeanMotion();
    i_xmdot_ = RecoveredMeanMotion() + 0.5 * temp1 * betao *
            i_x3thm1_ + 0.0625 * temp2 * betao *
            (13.0 - 78.0 * theta2 + 137.0 * theta4);
    double x1m5th = 1.0 - 5.0 * theta2;
    i_omgdot_ = -0.5 * temp1 * x1m5th +
            0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) +
            temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
    double xhdot1_ = -temp1 * i_cosio_;
    i_xnodot_ = xhdot1_ + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 *
            (3.0 - 7.0 * theta2)) * i_cosio_;
    i_xnodcf_ = 3.5 * betao2 * xhdot1_ * i_c1_;
    i_t2cof_ = 1.5 * i_c1_;

    if (fabs(i_cosio_ + 1.0) > 1.5e-12)
        i_xlcof_ = 0.125 * i_a3ovk2_ * i_sinio_ * (3.0 + 5.0 * i_cosio_) / (1.0 + i_cosio_);
    else
        i_xlcof_ = 0.125 * i_a3ovk2_ * i_sinio_ * (3.0 + 5.0 * i_cosio_) / 1.5e-12;

    i_aycof_ = 0.25 * i_a3ovk2_ * i_sinio_;
    i_x7thm1_ = 7.0 * theta2 - 1.0;

    if (i_use_deep_space_) {

        const double sing = sin(ArgumentPerigee());
        const double cosg = cos(ArgumentPerigee());
        i_gsto_ = Epoch().ToGMST();

        DeepSpaceInitialize(eosq, i_sinio_, i_cosio_, betao,
                theta2, sing, cosg, betao2,
                i_xmdot_, i_omgdot_, i_xnodot_);

    } else {

        double c3 = 0.0;
        if (Eccentricity() > 1.0e-4) {
            c3 = coef * tsi * i_a3ovk2_ * RecoveredMeanMotion() * constants_.AE *
                    i_sinio_ / Eccentricity();
        }

        i_c5_ = 2.0 * coef1 * RecoveredSemiMajorAxis() * betao2 * (1.0 + 2.75 *
                (etasq + eeta) + eeta * etasq);
        i_omgcof_ = BStar() * c3 * cos(ArgumentPerigee());

        i_xmcof_ = 0.0;
        if (Eccentricity() > 1.0e-4)
            i_xmcof_ = -constants_.TWOTHRD * coef * BStar() * constants_.AE / eeta;

        i_delmo_ = pow(1.0 + i_eta_ * (cos(MeanAnomoly())), 3.0);
        i_sinmo_ = sin(MeanAnomoly());

        if (!i_use_simple_model_) {
            double c1sq = i_c1_ * i_c1_;
            i_d2_ = 4.0 * RecoveredSemiMajorAxis() * tsi * c1sq;
            double temp = i_d2_ * tsi * i_c1_ / 3.0;
            i_d3_ = (17.0 * RecoveredSemiMajorAxis() + s4) * temp;
            i_d4_ = 0.5 * temp * RecoveredSemiMajorAxis() *
                    tsi * (221.0 * RecoveredSemiMajorAxis() + 31.0 * s4) * i_c1_;
            i_t3cof_ = i_d2_ + 2.0 * c1sq;
            i_t4cof_ = 0.25 * (3.0 * i_d3_ + i_c1_ *
                    (12.0 * i_d2_ + 10.0 * c1sq));
            i_t5cof_ = 0.2 * (3.0 * i_d4_ + 12.0 * i_c1_ *
                    i_d3_ + 6.0 * i_d2_ * i_d2_ + 15.0 *
                    c1sq * (2.0 * i_d2_ + c1sq));
        }
    }

    Eci eci;
    FindPosition(eci, 0.0);

    first_run_ = false;
}

void SGDP4::FindPosition(Eci& eci, double tsince) {

    if (i_use_deep_space_)
        FindPositionSDP4(eci, tsince);
    else
        FindPositionSGP4(eci, tsince);
}

void SGDP4::FindPositionSDP4(Eci& eci, double tsince) {

    /*
     * the final values
     */
    double e = Eccentricity();
    double a;
    double omega;
    double xl;
    double xnode;
    double xincl = Inclination();

    /*
     * update for secular gravity and atmospheric drag
     */
    double xmdf = MeanAnomoly() + i_xmdot_ * tsince;
    double omgadf = ArgumentPerigee() + i_omgdot_ * tsince;
    const double xnoddf = AscendingNode() + i_xnodot_ * tsince;

    omega = omgadf;

    const double tsq = tsince * tsince;
    xnode = xnoddf + i_xnodcf_ * tsq;
    double tempa = 1.0 - i_c1_ * tsince;
    double tempe = BStar() * i_c4_ * tsince;
    double templ = i_t2cof_ * tsq;

    double xn = RecoveredMeanMotion();

    /* t, xll, omgasm, xnodes, em, xinc, xn */
    DeepSpaceSecular(tsince, xmdf, omega, xnode, e, xincl, xn);

    if (xn <= 0.0) {
        throw new SatelliteException("Error: 2 (xn <= 0.0)");
    }

    a = pow(constants_.XKE / xn, constants_.TWOTHRD) * pow(tempa, 2.0);

    e -= tempe;
    /*
     * fix tolerance for error recognition
     */
    if (e >= 1.0 || e < -0.001) {
        throw new SatelliteException("Error: 1 (e >= 1.0 || e < -0.001)");
    }
    /*
     * fix tolerance to avoid a divide by zero
     */
    if (e < 1.0e-6)
        e = 1.0e-6;

    xmdf += RecoveredMeanMotion() * templ;
    double xlm = xmdf + omega + xnode;
    xnode = fmod(xnode, Globals::TWOPI());
    omega = fmod(omega, Globals::TWOPI());
    xlm = fmod(xlm, Globals::TWOPI());
    double xmam = fmod(xlm - omega - xnode, Globals::TWOPI());

    DeepSpacePeriodics(tsince, e, xincl, omega, xnode, xmam);

    /*
     * keeping xincl positive important unless you need to display xincl
     * and dislike negative inclinations
     */
    if (xincl < 0.0) {
        xincl = -xincl;
        xnode += Globals::PI();
        omega -= Globals::PI();
    }

    xl = xmam + omega + xnode;

    if (e < 0.0 || e > 1.0) {
        throw new SatelliteException("Error: 3 (e < 0.0 || e > 1.0)");
    }

    /*
     * re-compute the perturbed values
     */
    const double perturbed_sinio = sin(xincl);
    const double perturbed_cosio = cos(xincl);

    const double perturbed_theta2 = perturbed_cosio * perturbed_cosio;

    const double perturbed_x3thm1 = 3.0 * perturbed_theta2 - 1.0;
    const double perturbed_x1mth2 = 1.0 - perturbed_theta2;
    const double perturbed_x7thm1 = 7.0 * perturbed_theta2 - 1.0;

    double perturbed_xlcof;
    if (fabs(perturbed_cosio + 1.0) > 1.5e-12)
        perturbed_xlcof = 0.125 * i_a3ovk2_ * perturbed_sinio * (3.0 + 5.0 * perturbed_cosio) / (1.0 + perturbed_cosio);
    else
        perturbed_xlcof = 0.125 * i_a3ovk2_ * perturbed_sinio * (3.0 + 5.0 * perturbed_cosio) / 1.5e-12;

    const double perturbed_aycof = 0.25 * i_a3ovk2_ * perturbed_sinio;

    /*
     * using calculated values, find position and velocity
     */
    CalculateFinalPositionVelocity(eci, tsince, e,
            a, omega, xl, xnode,
            xincl, perturbed_xlcof, perturbed_aycof,
            perturbed_x3thm1, perturbed_x1mth2, perturbed_x7thm1,
            perturbed_cosio, perturbed_sinio);

}

void SGDP4::FindPositionSGP4(Eci& eci, double tsince) const {

    /*
     * the final values
     */
    double e;
    double a;
    double omega;
    double xl;
    double xnode;
    double xincl;

    /*
     * update for secular gravity and atmospheric drag
     */
    const double xmdf = MeanAnomoly() + i_xmdot_ * tsince;
    const double omgadf = ArgumentPerigee() + i_omgdot_ * tsince;
    const double xnoddf = AscendingNode() + i_xnodot_ * tsince;

    const double tsq = tsince * tsince;
    xnode = xnoddf + i_xnodcf_ * tsq;
    double tempa = 1.0 - i_c1_ * tsince;
    double tempe = BStar() * i_c4_ * tsince;
    double templ = i_t2cof_ * tsq;

    xincl = Inclination();
    omega = omgadf;
    double xmp = xmdf;

    if (!i_use_simple_model_) {
        const double delomg = i_omgcof_ * tsince;
        const double delm = i_xmcof_ * (pow(1.0 + i_eta_ * cos(xmdf), 3.0) - i_delmo_);
        const double temp = delomg + delm;

        xmp += temp;
        omega -= temp;

        const double tcube = tsq * tsince;
        const double tfour = tsince * tcube;

        tempa -= i_d2_ * tsq - i_d3_ * tcube - i_d4_ * tfour;
        tempe += BStar() * i_c5_ * (sin(xmp) - i_sinmo_);
        templ += i_t3cof_ * tcube + tfour * (i_t4cof_ + tsince * i_t5cof_);
    }

    a = RecoveredSemiMajorAxis() * pow(tempa, 2.0);
    e = Eccentricity() - tempe;
    xl = xmp + omega + xnode + RecoveredMeanMotion() * templ;

    /*
     * using calculated values, find position and velocity
     * we can pass in constants from Initialize() as these dont change
     */
    CalculateFinalPositionVelocity(eci, tsince, e,
            a, omega, xl, xnode,
            xincl, i_xlcof_, i_aycof_,
            i_x3thm1_, i_x1mth2_, i_x7thm1_,
            i_cosio_, i_sinio_);

}

void SGDP4::CalculateFinalPositionVelocity(Eci& eci, const double& tsince, const double& e,
        const double& a, const double& omega, const double& xl, const double& xnode,
        const double& xincl, const double& xlcof, const double& aycof,
        const double& x3thm1, const double& x1mth2, const double& x7thm1,
        const double& cosio, const double& sinio) const {

    double temp;
    double temp1;
    double temp2;
    double temp3;

    if (a < 1.0) {
        throw new SatelliteException("Error: Satellite crashed (a < 1.0)");
    }

    if (e < -1.0e-3) {
        throw new SatelliteException("Error: Modified eccentricity too low (e < -1.0e-3)");
    }

    const double beta = sqrt(1.0 - e * e);
    const double xn = constants_.XKE / pow(a, 1.5);
    /*
     * long period periodics
     */
    const double axn = e * cos(omega);
    temp = 1.0 / (a * beta * beta);
    const double xll = temp * xlcof * axn;
    const double aynl = temp * aycof;
    const double xlt = xl + xll;
    const double ayn = e * sin(omega) + aynl;
    const double elsq = axn * axn + ayn * ayn;

    if (elsq >= 1.0) {
        throw new SatelliteException("Error: sqrt(e) >= 1 (elsq >= 1.0)");
    }

    /*
     * solve keplers equation
     * - solve using Newton-Raphson root solving
     * - here capu is almost the mean anomoly
     * - initialise the eccentric anomaly term epw
     * - The fmod saves reduction of angle to +/-2pi in sin/cos() and prevents
     * convergence problems.
     */
    const double capu = fmod(xlt - xnode, Globals::TWOPI());
    double epw = capu;

    double sinepw = 0.0;
    double cosepw = 0.0;
    double ecose = 0.0;
    double esine = 0.0;

    /*
     * sensibility check for N-R correction
     */
    const double max_newton_naphson = 1.25 * fabs(sqrt(elsq));

    bool kepler_running = true;

    for (int i = 0; i < 10 && kepler_running; i++) {
        sinepw = sin(epw);
        cosepw = cos(epw);
        ecose = axn * cosepw + ayn * sinepw;
        esine = axn * sinepw - ayn * cosepw;

        double f = capu - epw + esine;

        if (fabs(f) < 1.0e-12) {
            kepler_running = false;
        } else {
            /*
             * 1st order Newton-Raphson correction
             */
            const double fdot = 1.0 - ecose;
            double delta_epw = f / fdot;

            /*
             * 2nd order Newton-Raphson correction.
             * f / (fdot - 0.5 * d2f * f/fdot)
             */
            if (i == 0) {
                if (delta_epw > max_newton_naphson)
                    delta_epw = max_newton_naphson;
                else if (delta_epw < -max_newton_naphson)
                    delta_epw = -max_newton_naphson;
            } else {
                delta_epw = f / (fdot + 0.5 * esine * delta_epw);
            }

            /*
             * Newton-Raphson correction of -F/DF
             */
            epw += delta_epw;
        }
    }
    /*
     * short period preliminary quantities
     */
    temp = 1.0 - elsq;
    const double pl = a * temp;
    const double r = a * (1.0 - ecose);
    temp1 = 1.0 / r;
    const double rdot = constants_.XKE * sqrt(a) * esine * temp1;
    const double rfdot = constants_.XKE * sqrt(pl) * temp1;
    temp2 = a * temp1;
    const double betal = sqrt(temp);
    temp3 = 1.0 / (1.0 + betal);
    const double cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
    const double sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
    const double u = atan2(sinu, cosu);
    const double sin2u = 2.0 * sinu * cosu;
    const double cos2u = 2.0 * cosu * cosu - 1.0;
    temp = 1.0 / pl;
    temp1 = constants_.CK2 * temp;
    temp2 = temp1 * temp;

    /*
     * update for short periodics
     */
    const double rk = r * (1.0 - 1.5 * temp2 * betal * x3thm1) + 0.5 * temp1 * x1mth2 * cos2u;
    const double uk = u - 0.25 * temp2 * x7thm1 * sin2u;
    const double xnodek = xnode + 1.5 * temp2 * cosio * sin2u;
    const double xinck = xincl + 1.5 * temp2 * cosio * sinio * cos2u;
    const double rdotk = rdot - xn * temp1 * x1mth2 * sin2u;
    const double rfdotk = rfdot + xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1);

    if (rk < 0.0) {
        throw new SatelliteException("Error: satellite decayed (rk < 0.0)");
    }

    /*
     * orientation vectors
     */
    const double sinuk = sin(uk);
    const double cosuk = cos(uk);
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
    const double x = rk * ux * constants_.XKMPER;
    const double y = rk * uy * constants_.XKMPER;
    const double z = rk * uz * constants_.XKMPER;
    Vector position(x, y, z);
    const double xdot = (rdotk * ux + rfdotk * vx) * constants_.XKMPER / 60.0;
    const double ydot = (rdotk * uy + rfdotk * vy) * constants_.XKMPER / 60.0;
    const double zdot = (rdotk * uz + rfdotk * vz) * constants_.XKMPER / 60.0;
    Vector velocity(xdot, ydot, zdot);

    Julian julian = Epoch();
    julian.AddMin(tsince);
    eci = Eci(julian, position, velocity);
}

/*
 * deep space initialization
 */
void SGDP4::DeepSpaceInitialize(const double& eosq, const double& sinio, const double& cosio, const double& betao,
        const double& theta2, const double& sing, const double& cosg, const double& betao2,
        const double& xmdot, const double& omgdot, const double& xnodot) {

    double se = 0.0;
    double si = 0.0;
    double sl = 0.0;
    double sgh = 0.0;
    double shdq = 0.0;

    double bfact = 0.0;

    static const double ZNS = 1.19459E-5;
    static const double C1SS = 2.9864797E-6;
    static const double ZES = 0.01675;
    static const double ZNL = 1.5835218E-4;
    static const double C1L = 4.7968065E-7;
    static const double ZEL = 0.05490;
    static const double ZCOSIS = 0.91744867;
    static const double ZSINI = 0.39785416;
    static const double ZSINGS = -0.98088458;
    static const double ZCOSGS = 0.1945905;
    static const double Q22 = 1.7891679E-6;
    static const double Q31 = 2.1460748E-6;
    static const double Q33 = 2.2123015E-7;
    static const double ROOT22 = 1.7891679E-6;
    static const double ROOT32 = 3.7393792E-7;
    static const double ROOT44 = 7.3636953E-9;
    static const double ROOT52 = 1.1428639E-7;
    static const double ROOT54 = 2.1765803E-9;
    static const double THDT = 4.3752691E-3;

    const double aqnv = 1.0 / RecoveredSemiMajorAxis();
    const double xpidot = omgdot + xnodot;
    const double sinq = sin(AscendingNode());
    const double cosq = cos(AscendingNode());

    /*
     * initialize lunar / solar terms
     */
    const double d_day_ = Epoch().FromJan1_12h_1900();

    const double xnodce = fmod(4.5236020 - 9.2422029e-4 * d_day_, Globals::TWOPI());
    const double stem = sin(xnodce);
    const double ctem = cos(xnodce);
    const double zcosil = 0.91375164 - 0.03568096 * ctem;
    const double zsinil = sqrt(1.0 - zcosil * zcosil);
    const double zsinhl = 0.089683511 * stem / zsinil;
    const double zcoshl = sqrt(1.0 - zsinhl * zsinhl);
    const double c = 4.7199672 + 0.22997150 * d_day_;
    const double gam = 5.8351514 + 0.0019443680 * d_day_;
    d_zmol_ = Globals::Fmod2p(c - gam);
    double zx = 0.39785416 * stem / zsinil;
    double zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    /*
     * todo: check
     */
    zx = atan2(zx, zy);
    zx = fmod(gam + zx - xnodce, Globals::TWOPI());

    const double zcosgl = cos(zx);
    const double zsingl = sin(zx);
    d_zmos_ = 6.2565837 + 0.017201977 * d_day_;
    d_zmos_ = Globals::Fmod2p(d_zmos_);

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
    const double xnoi = 1.0 / RecoveredMeanMotion();

    for (int cnt = 0; cnt < 2; cnt++) {
        /*
         * solar terms are done a second time after lunar terms are done
         */
        const double a1 = zcosg * zcosh + zsing * zcosi * zsinh;
        const double a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
        const double a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
        const double a8 = zsing * zsini;
        const double a9 = zsing * zsinh + zcosg * zcosi*zcosh;
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
        const double s1 = -15.0 * Eccentricity() * s4;
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
        if (Inclination() < 5.2359877e-2 || Inclination() > Globals::PI() - 5.2359877e-2) {
            shdq = 0.0;
        } else {
            shdq = (-zn * s2 * (z21 + z23)) / sinio;
        }

        d_ee2_ = 2.0 * s1 * s6;
        d_e3_ = 2.0 * s1 * s7;
        d_xi2_ = 2.0 * s2 * z12;
        d_xi3_ = 2.0 * s2 * (z13 - z11);
        d_xl2_ = -2.0 * s3 * z2;
        d_xl3_ = -2.0 * s3 * (z3 - z1);
        d_xl4_ = -2.0 * s3 * (-21.0 - 9.0 * eosq) * ze;
        d_xgh2_ = 2.0 * s4 * z32;
        d_xgh3_ = 2.0 * s4 * (z33 - z31);
        d_xgh4_ = -18.0 * s4 * ze;
        d_xh2_ = -2.0 * s2 * z22;
        d_xh3_ = -2.0 * s2 * (z23 - z21);

        if (cnt == 1)
            break;
        /*
         * do lunar terms
         */
        d_sse_ = se;
        d_ssi_ = si;
        d_ssl_ = sl;
        d_ssh_ = shdq;
        d_ssg_ = sgh - cosio * d_ssh_;
        d_se2_ = d_ee2_;
        d_si2_ = d_xi2_;
        d_sl2_ = d_xl2_;
        d_sgh2_ = d_xgh2_;
        d_sh2_ = d_xh2_;
        d_se3_ = d_e3_;
        d_si3_ = d_xi3_;
        d_sl3_ = d_xl3_;
        d_sgh3_ = d_xgh3_;
        d_sh3_ = d_xh3_;
        d_sl4_ = d_xl4_;
        d_sgh4_ = d_xgh4_;
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

    d_sse_ += se;
    d_ssi_ += si;
    d_ssl_ += sl;
    d_ssg_ += sgh - cosio * shdq;
    d_ssh_ += shdq;


    d_resonance_flag_ = false;
    d_synchronous_flag_ = false;
    bool initialize_integrator = true;

    if (RecoveredMeanMotion() < 0.0052359877 && RecoveredMeanMotion() > 0.0034906585) {

        /*
         * 24h synchronous resonance terms initialization
         */
        d_resonance_flag_ = true;
        d_synchronous_flag_ = true;

        const double g200 = 1.0 + eosq * (-2.5 + 0.8125 * eosq);
        const double g310 = 1.0 + 2.0 * eosq;
        const double g300 = 1.0 + eosq * (-6.0 + 6.60937 * eosq);
        const double f220 = 0.75 * (1.0 + cosio) * (1.0 + cosio);
        const double f311 = 0.9375 * sinio * sinio * (1.0 + 3.0 * cosio) - 0.75 * (1.0 + cosio);
        double f330 = 1.0 + cosio;
        f330 = 1.875 * f330 * f330 * f330;
        d_del1_ = 3.0 * RecoveredMeanMotion() * RecoveredMeanMotion() * aqnv * aqnv;
        d_del2_ = 2.0 * d_del1_ * f220 * g200 * Q22;
        d_del3_ = 3.0 * d_del1_ * f330 * g300 * Q33 * aqnv;
        d_del1_ = d_del1_ * f311 * g310 * Q31 * aqnv;
        d_fasx2_ = 0.13130908;
        d_fasx4_ = 2.8843198;
        d_fasx6_ = 0.37448087;

        d_xlamo_ = MeanAnomoly() + AscendingNode() + ArgumentPerigee() - i_gsto_;
        bfact = xmdot + xpidot - THDT;
        bfact += d_ssl_ + d_ssg_ + d_ssh_;

    } else if (RecoveredMeanMotion() < 8.26e-3 || RecoveredMeanMotion() > 9.24e-3 || Eccentricity() < 0.5) {
        initialize_integrator = false;
    } else {
        /*
         * geopotential resonance initialization for 12 hour orbits
         */
        d_resonance_flag_ = true;

        const double eoc = Eccentricity() * eosq;

        double g211;
        double g310;
        double g322;
        double g410;
        double g422;
        double g520;

        double g201 = -0.306 - (Eccentricity() - 0.64) * 0.440;

        if (Eccentricity() <= 0.65) {
            g211 = 3.616 - 13.247 * Eccentricity() + 16.290 * eosq;
            g310 = -19.302 + 117.390 * Eccentricity() - 228.419 * eosq + 156.591 * eoc;
            g322 = -18.9068 + 109.7927 * Eccentricity() - 214.6334 * eosq + 146.5816 * eoc;
            g410 = -41.122 + 242.694 * Eccentricity() - 471.094 * eosq + 313.953 * eoc;
            g422 = -146.407 + 841.880 * Eccentricity() - 1629.014 * eosq + 1083.435 * eoc;
            g520 = -532.114 + 3017.977 * Eccentricity() - 5740 * eosq + 3708.276 * eoc;
        } else {
            g211 = -72.099 + 331.819 * Eccentricity() - 508.738 * eosq + 266.724 * eoc;
            g310 = -346.844 + 1582.851 * Eccentricity() - 2415.925 * eosq + 1246.113 * eoc;
            g322 = -342.585 + 1554.908 * Eccentricity() - 2366.899 * eosq + 1215.972 * eoc;
            g410 = -1052.797 + 4758.686 * Eccentricity() - 7193.992 * eosq + 3651.957 * eoc;
            g422 = -3581.69 + 16178.11 * Eccentricity() - 24462.77 * eosq + 12422.52 * eoc;

            if (Eccentricity() <= 0.715) {
                g520 = 1464.74 - 4664.75 * Eccentricity() + 3763.64 * eosq;
            } else {
                g520 = -5149.66 + 29936.92 * Eccentricity() - 54087.36 * eosq + 31324.56 * eoc;
            }
        }

        double g533;
        double g521;
        double g532;

        if (Eccentricity() < 0.7) {
            g533 = -919.2277 + 4988.61 * Eccentricity() - 9064.77 * eosq + 5542.21 * eoc;
            g521 = -822.71072 + 4568.6173 * Eccentricity() - 8491.4146 * eosq + 5337.524 * eoc;
            g532 = -853.666 + 4690.25 * Eccentricity() - 8624.77 * eosq + 5341.4 * eoc;
        } else {
            g533 = -37995.78 + 161616.52 * Eccentricity() - 229838.2 * eosq + 109377.94 * eoc;
            g521 = -51752.104 + 218913.95 * Eccentricity() - 309468.16 * eosq + 146349.42 * eoc;
            g532 = -40023.88 + 170470.89 * Eccentricity() - 242699.48 * eosq + 115605.82 * eoc;
        }

        const double sini2 = sinio * sinio;
        const double f220 = 0.75 * (1.0 + 2.0 * cosio + theta2);
        const double f221 = 1.5 * sini2;
        const double f321 = 1.875 * sinio * (1.0 - 2.0 * cosio - 3.0 * theta2);
        const double f322 = -1.875 * sinio * (1.0 + 2.0 * cosio - 3.0 * theta2);
        const double f441 = 35.0 * sini2 * f220;
        const double f442 = 39.3750 * sini2 * sini2;
        const double f522 = 9.84375 * sinio * (sini2 * (1.0 - 2.0 * cosio - 5.0 * theta2)
                + 0.33333333 * (-2.0 + 4.0 * cosio + 6.0 * theta2));
        const double f523 = sinio * (4.92187512 * sini2 * (-2.0 - 4.0 * cosio + 10.0 * theta2)
                + 6.56250012 * (1.0 + 2.0 * cosio - 3.0 * theta2));
        const double f542 = 29.53125 * sinio * (2.0 - 8.0 * cosio + theta2 *
                (-12.0 + 8.0 * cosio + 10.0 * theta2));
        const double f543 = 29.53125 * sinio * (-2.0 - 8.0 * cosio + theta2 *
                (12.0 + 8.0 * cosio - 10.0 * theta2));

        const double xno2 = RecoveredMeanMotion() * RecoveredMeanMotion();
        const double ainv2 = aqnv * aqnv;

        double temp1 = 3.0 * xno2 * ainv2;
        double temp = temp1 * ROOT22;
        d_d2201_ = temp * f220 * g201;
        d_d2211_ = temp * f221 * g211;
        temp1 = temp1 * aqnv;
        temp = temp1 * ROOT32;
        d_d3210_ = temp * f321 * g310;
        d_d3222_ = temp * f322 * g322;
        temp1 = temp1 * aqnv;
        temp = 2.0 * temp1 * ROOT44;
        d_d4410_ = temp * f441 * g410;
        d_d4422_ = temp * f442 * g422;
        temp1 = temp1 * aqnv;
        temp = temp1 * ROOT52;
        d_d5220_ = temp * f522 * g520;
        d_d5232_ = temp * f523 * g532;
        temp = 2.0 * temp1 * ROOT54;
        d_d5421_ = temp * f542 * g521;
        d_d5433_ = temp * f543 * g533;

        d_xlamo_ = MeanAnomoly() + AscendingNode() + AscendingNode() - i_gsto_ - i_gsto_;
        bfact = xmdot + xnodot + xnodot - THDT - THDT;
        bfact = bfact + d_ssl_ + d_ssh_ + d_ssh_;
    }

    if (initialize_integrator) {
        /*
         * initialize integrator
         */
        d_xfact_ = bfact - RecoveredMeanMotion();
        d_atime_ = 0.0;
        d_xni_ = RecoveredMeanMotion();
        d_xli_ = d_xlamo_;
    }
}

void SGDP4::DeepSpaceCalculateLunarSolarTerms(const double t, double& pe, double& pinc,
        double& pl, double& pgh, double& ph) const {

    static const double ZES = 0.01675;
    static const double ZNS = 1.19459E-5;
    static const double ZNL = 1.5835218E-4;
    static const double ZEL = 0.05490;

    /*
     * calculate solar terms
     */
    double zm = d_zmos_ + ZNS * t;
    if (first_run_)
        zm = d_zmos_;
    double zf = zm + 2.0 * ZES * sin(zm);
    double sinzf = sin(zf);
    double f2 = 0.5 * sinzf * sinzf - 0.25;
    double f3 = -0.5 * sinzf * cos(zf);
    const double ses = d_se2_ * f2 + d_se3_ * f3;
    const double sis = d_si2_ * f2 + d_si3_ * f3;
    const double sls = d_sl2_ * f2 + d_sl3_ * f3 + d_sl4_ * sinzf;
    const double sghs = d_sgh2_ * f2 + d_sgh3_ * f3 + d_sgh4_ * sinzf;
    const double shs = d_sh2_ * f2 + d_sh3_ * f3;

    /*
     * calculate lunar terms
     */
    zm = d_zmol_ + ZNL * t;
    if (first_run_)
        zm = d_zmol_;
    zf = zm + 2.0 * ZEL * sin(zm);
    sinzf = sin(zf);
    f2 = 0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * cos(zf);
    const double sel = d_ee2_ * f2 + d_e3_ * f3;
    const double sil = d_xi2_ * f2 + d_xi3_ * f3;
    const double sll = d_xl2_ * f2 + d_xl3_ * f3 + d_xl4_ * sinzf;
    const double sghl = d_xgh2_ * f2 + d_xgh3_ * f3 + d_xgh4_ * sinzf;
    const double shl = d_xh2_ * f2 + d_xh3_ * f3;

    /*
     * merge calculated values
     */
    pe = ses + sel;
    pinc = sis + sil;
    pl = sls + sll;
    pgh = sghs + sghl;
    ph = shs + shl;
}

/*
 * calculate lunar / solar periodics and apply
 */
void SGDP4::DeepSpacePeriodics(const double& t, double& em,
        double& xinc, double& omgasm, double& xnodes, double& xll) {

    /*
     * storage for lunar / solar terms set by DeepSpaceCalculateLunarSolarTerms()
     */
    double pe = 0.0;
    double pinc = 0.0;
    double pl = 0.0;
    double pgh = 0.0;
    double ph = 0.0;

    /*
     * calculate lunar / solar terms for current time
     */
    DeepSpaceCalculateLunarSolarTerms(t, pe, pinc, pl, pgh, ph);

    if (!first_run_) {

        xinc += pinc;
        em += pe;

        /* Spacetrack report #3 has sin/cos from before perturbations
         * added to xinc (oldxinc), but apparently report # 6 has then
         * from after they are added.
         * use for strn3
         * if (Inclination() >= 0.2)
         * use for gsfc
         * if (xinc >= 0.2)
         * (moved from start of function)
         */
        const double sinis = sin(xinc);
        const double cosis = cos(xinc);

        if (xinc >= 0.2) {
            /*
             * apply periodics directly
             */
            const double tmp_ph = ph / sinis;

            omgasm += pgh - cosis * tmp_ph;
            xnodes += tmp_ph;
            xll += pl;
        } else {
            /*
             * apply periodics with lyddane modification
             */
            const double sinok = sin(xnodes);
            const double cosok = cos(xnodes);
            double alfdp = sinis * sinok;
            double betdp = sinis * cosok;
            const double dalf = ph * cosok + pinc * cosis * sinok;
            const double dbet = -ph * sinok + pinc * cosis * cosok;

            alfdp += dalf;
            betdp += dbet;

            double xls = xll + omgasm + cosis * xnodes;
            double dls = pl + pgh - pinc * xnodes * sinis;
            xls += dls;

            /*
             * save old xnodes value
             */
            const double oldxnodes = xnodes;

            xnodes = atan2(alfdp, betdp);
            /*
             * Get perturbed xnodes in to same quadrant as original.
             */
            if (fabs(oldxnodes - xnodes) > Globals::PI()) {
                if (xnodes < oldxnodes)
                    xnodes += Globals::TWOPI();
                else
                    xnodes -= Globals::TWOPI();
            }

            xll += pl;
            omgasm = xls - xll - cosis * xnodes;
        }
    }
}

/*
 * deep space secular effects
 */
void SGDP4::DeepSpaceSecular(const double& t, double& xll, double& omgasm,
        double& xnodes, double& em, double& xinc, double& xn) const {

    static const double THDT = 4.37526908801129966e-3;

    static const double STEP2 = 259200.0;
    static const double STEP = 720.0;

    double xldot = 0.0;
    double xndot = 0.0;
    double xnddt = 0.0;

    xll += d_ssl_ * t;
    omgasm += d_ssg_ * t;
    xnodes += d_ssh_ * t;
    em += d_sse_ * t;
    xinc += d_ssi_ * t;

    if (!d_resonance_flag_)
        return;

    /*
     * 1st condition (if t is less than one time step from epoch)
     * 2nd condition (if d_atime_ and t are of opposite signs, so zero crossing required)
     * 3rd condition (if t is closer to zero than d_atime_)
     */
    if (fabs(t) < STEP ||
            t * d_atime_ <= 0.0 ||
            fabs(t) < fabs(d_atime_)) {
        /*
         * restart from epoch
         */
        d_atime_ = 0.0;
        d_xni_ = RecoveredMeanMotion();
        d_xli_ = d_xlamo_;
    }

    double ft = t - d_atime_;

    /*
     * if time difference (ft) is greater than the time step (720.0)
     * loop around until d_atime_ is within one time step of t
     */
    if (fabs(ft) >= STEP) {

        /*
         * calculate step direction to allow d_atime_
         * to catch up with t
         */
        double delt = -STEP;
        if (ft >= 0.0)
            delt = STEP;

        do {
            DeepSpaceCalcIntegrator(delt, STEP2, xndot, xnddt, xldot);

            ft = t - d_atime_;
        } while (fabs(ft) >= STEP);
    }

    /*
     * calculate dot terms
     */
    DeepSpaceCalcDotTerms(xndot, xnddt, xldot);

    /*
     * integrator
     */
    xn = d_xni_ + xndot * ft + xnddt * ft * ft * 0.5;
    const double xl = d_xli_ + xldot * ft + xndot * ft * ft * 0.5;
    const double temp = -xnodes + i_gsto_ + t * THDT;

    if (d_synchronous_flag_)
        xll = xl + temp - omgasm;
    else
        xll = xl + temp + temp;
}

/*
 * calculate dot terms
 */
void SGDP4::DeepSpaceCalcDotTerms(double& xndot, double& xnddt, double& xldot) const {

    static const double G22 = 5.7686396;
    static const double G32 = 0.95240898;
    static const double G44 = 1.8014998;
    static const double G52 = 1.0508330;
    static const double G54 = 4.4108898;

    if (d_synchronous_flag_) {

        xndot = d_del1_ * sin(d_xli_ - d_fasx2_) +
                d_del2_ * sin(2.0 * (d_xli_ - d_fasx4_)) +
                d_del3_ * sin(3.0 * (d_xli_ - d_fasx6_));
        xnddt = d_del1_ * cos(d_xli_ - d_fasx2_) + 2.0 *
                d_del2_ * cos(2.0 * (d_xli_ - d_fasx4_)) + 3.0 *
                d_del3_ * cos(3.0 * (d_xli_ - d_fasx6_));

    } else {

        const double xomi = ArgumentPerigee() + i_omgdot_ * d_atime_;
        const double x2omi = xomi + xomi;
        const double x2li = d_xli_ + d_xli_;

        xndot = d_d2201_ * sin(x2omi + d_xli_ - G22)
                + d_d2211_ * sin(d_xli_ - G22)
                + d_d3210_ * sin(xomi + d_xli_ - G32)
                + d_d3222_ * sin(-xomi + d_xli_ - G32)
                + d_d4410_ * sin(x2omi + x2li - G44)
                + d_d4422_ * sin(x2li - G44)
                + d_d5220_ * sin(xomi + d_xli_ - G52)
                + d_d5232_ * sin(-xomi + d_xli_ - G52)
                + d_d5421_ * sin(xomi + x2li - G54)
                + d_d5433_ * sin(-xomi + x2li - G54);
        xnddt = d_d2201_ * cos(x2omi + d_xli_ - G22)
                + d_d2211_ * cos(d_xli_ - G22)
                + d_d3210_ * cos(xomi + d_xli_ - G32)
                + d_d3222_ * cos(-xomi + d_xli_ - G32)
                + d_d5220_ * cos(xomi + d_xli_ - G52)
                + d_d5232_ * cos(-xomi + d_xli_ - G52)
                + 2.0 * (d_d4410_ * cos(x2omi + x2li - G44)
                + d_d4422_ * cos(x2li - G44)
                + d_d5421_ * cos(xomi + x2li - G54)
                + d_d5433_ * cos(-xomi + x2li - G54));
    }

    xldot = d_xni_ + d_xfact_;
    xnddt *= xldot;
}

/*
 * deep space integrator for time period of delt
 */
void SGDP4::DeepSpaceCalcIntegrator(const double& delt, const double& step2, double& xndot, double& xnddt, double& xldot) const {

    /*
     * calculate dot terms
     */
    DeepSpaceCalcDotTerms(xndot, xnddt, xldot);

    /*
     * integrator
     */
    d_xli_ += xldot * delt + xndot * step2;
    d_xni_ += xndot * delt + xnddt * step2;

    /*
     * increment integrator time
     */
    d_atime_ += delt;
}
