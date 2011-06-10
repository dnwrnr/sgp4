#include "SGP4.h"

#include "Vector.h"
#include "SatelliteException.h"

#include <cmath>
#include <iomanip>
#include <cstring>

SGP4::SGP4(void) {

    Reset();
}

SGP4::SGP4(const Tle& tle) {

    SetTle(tle);
}

SGP4::~SGP4(void) {
}

void SGP4::SetTle(const Tle& tle) {

    /*
     * reset all constants etc
     */
    Reset();

    /*
     * extract and format tle data
     */
    mean_anomoly_ = tle.MeanAnomaly(false);
    ascending_node_ = tle.RightAscendingNode(false);
    argument_perigee_ = tle.ArgumentPerigee(false);
    eccentricity_ = tle.Eccentricity();
    inclination_ = tle.Inclination(false);
    mean_motion_ = tle.MeanMotion() * kTWOPI / kMINUTES_PER_DAY;
    bstar_ = tle.BStar();
    epoch_ = tle.Epoch();

    /*
     * error checks
     */
    if (eccentricity_ < 0.0 || eccentricity_ > 1.0 - 1.0e-3) {
        throw SatelliteException("Eccentricity out of range");
    }

    if (inclination_ < 0.0 || eccentricity_ > kPI) {
        throw SatelliteException("Inclination out of range");
    }

    Initialize();
}

void SGP4::Initialize() {






    /*
     * recover original mean motion (xnodp) and semimajor axis (aodp)
     * from input elements
     */
    const double a1 = pow(kXKE / MeanMotion(), kTWOTHIRD);
    common_consts_.cosio = cos(Inclination());
    common_consts_.sinio = sin(Inclination());
    const double theta2 = common_consts_.cosio * common_consts_.cosio;
    common_consts_.x3thm1 = 3.0 * theta2 - 1.0;
    const double eosq = Eccentricity() * Eccentricity();
    const double betao2 = 1.0 - eosq;
    const double betao = sqrt(betao2);
    const double temp = (1.5 * kCK2) * common_consts_.x3thm1 / (betao * betao2);
    const double del1 = temp / (a1 * a1);
    const double a0 = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + del1 * 134.0 / 81.0)));
    const double del0 = temp / (a0 * a0);

    recovered_mean_motion_ = MeanMotion() / (1.0 + del0);
    /*
     * alternative way to calculate
     * doesnt affect final results
     * recovered_semi_major_axis_ = pow(XKE / RecoveredMeanMotion(), TWOTHIRD);
     */
    recovered_semi_major_axis_ = a0 / (1.0 - del0);

    /*
     * find perigee and period
     */
    perigee_ = (RecoveredSemiMajorAxis() * (1.0 - Eccentricity()) - kAE) * kXKMPER;
    period_ = kTWOPI / RecoveredMeanMotion();









    if (Period() >= 225.0) {
        use_deep_space_ = true;
    } else {
        use_deep_space_ = false;
        use_simple_model_ = false;
        /*
         * for perigee less than 220 kilometers, the simple_model flag is set and
         * the equations are truncated to linear variation in sqrt a and
         * quadratic variation in mean anomly. also, the c3 term, the
         * delta omega term and the delta m term are dropped
         */
        if (Perigee() < 220.0) {
            use_simple_model_ = true;
        }
    }

    /*
     * for perigee below 156km, the values of
     * s4 and qoms2t are altered
     */
    double s4 = kS;
    double qoms24 = kQOMS2T;
    if (Perigee() < 156.0) {
        s4 = Perigee() - 78.0;
        if (Perigee() < 98.0) {
            s4 = 20.0;
        }
        qoms24 = pow((120.0 - s4) * kAE / kXKMPER, 4.0);
        s4 = s4 / kXKMPER + kAE;
    }

    /*
     * generate constants
     */
    const double pinvsq = 1.0 / (RecoveredSemiMajorAxis() * RecoveredSemiMajorAxis() * betao2 * betao2);
    const double tsi = 1.0 / (RecoveredSemiMajorAxis() - s4);
    common_consts_.eta = RecoveredSemiMajorAxis() * Eccentricity() * tsi;
    const double etasq = common_consts_.eta * common_consts_.eta;
    const double eeta = Eccentricity() * common_consts_.eta;
    const double psisq = fabs(1.0 - etasq);
    const double coef = qoms24 * pow(tsi, 4.0);
    const double coef1 = coef / pow(psisq, 3.5);
    const double c2 = coef1 * RecoveredMeanMotion() * (RecoveredSemiMajorAxis() *
            (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
            0.75 * kCK2 * tsi / psisq *
            common_consts_.x3thm1 * (8.0 + 3.0 * etasq *
            (8.0 + etasq)));
    common_consts_.c1 = BStar() * c2;
    common_consts_.a3ovk2 = -kXJ3 / kCK2 * pow(kAE, 3.0);
    common_consts_.x1mth2 = 1.0 - theta2;
    common_consts_.c4 = 2.0 * RecoveredMeanMotion() * coef1 * RecoveredSemiMajorAxis() * betao2 *
            (common_consts_.eta * (2.0 + 0.5 * etasq) + Eccentricity() * (0.5 + 2.0 * etasq) -
            2.0 * kCK2 * tsi / (RecoveredSemiMajorAxis() * psisq) *
            (-3.0 * common_consts_.x3thm1 * (1.0 - 2.0 * eeta + etasq *
            (1.5 - 0.5 * eeta)) + 0.75 * common_consts_.x1mth2 * (2.0 * etasq - eeta *
            (1.0 + etasq)) * cos(2.0 * ArgumentPerigee())));
    const double theta4 = theta2 * theta2;
    const double temp1 = 3.0 * kCK2 * pinvsq * RecoveredMeanMotion();
    const double temp2 = temp1 * kCK2 * pinvsq;
    const double temp3 = 1.25 * kCK4 * pinvsq * pinvsq * RecoveredMeanMotion();
    common_consts_.xmdot = RecoveredMeanMotion() + 0.5 * temp1 * betao *
            common_consts_.x3thm1 + 0.0625 * temp2 * betao *
            (13.0 - 78.0 * theta2 + 137.0 * theta4);
    const double x1m5th = 1.0 - 5.0 * theta2;
    common_consts_.omgdot = -0.5 * temp1 * x1m5th +
            0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) +
            temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
    const double xhdot1 = -temp1 * common_consts_.cosio;
    common_consts_.xnodot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 *
            (3.0 - 7.0 * theta2)) * common_consts_.cosio;
    common_consts_.xnodcf = 3.5 * betao2 * xhdot1 * common_consts_.c1;
    common_consts_.t2cof = 1.5 * common_consts_.c1;

    if (fabs(common_consts_.cosio + 1.0) > 1.5e-12)
        common_consts_.xlcof = 0.125 * common_consts_.a3ovk2 * common_consts_.sinio * (3.0 + 5.0 * common_consts_.cosio) / (1.0 + common_consts_.cosio);
    else
        common_consts_.xlcof = 0.125 * common_consts_.a3ovk2 * common_consts_.sinio * (3.0 + 5.0 * common_consts_.cosio) / 1.5e-12;

    common_consts_.aycof = 0.25 * common_consts_.a3ovk2 * common_consts_.sinio;
    common_consts_.x7thm1 = 7.0 * theta2 - 1.0;

    if (use_deep_space_) {

        deepspace_consts_.gsto = Epoch().ToGreenwichSiderealTime();

        DeepSpaceInitialize(eosq, common_consts_.sinio, common_consts_.cosio, betao,
                theta2, betao2,
                common_consts_.xmdot, common_consts_.omgdot, common_consts_.xnodot);

    } else {

        double c3 = 0.0;
        if (Eccentricity() > 1.0e-4) {
            c3 = coef * tsi * common_consts_.a3ovk2 * RecoveredMeanMotion() * kAE *
                    common_consts_.sinio / Eccentricity();
        }

        nearspace_consts_.c5 = 2.0 * coef1 * RecoveredSemiMajorAxis() * betao2 * (1.0 + 2.75 *
                (etasq + eeta) + eeta * etasq);
        nearspace_consts_.omgcof = BStar() * c3 * cos(ArgumentPerigee());

        nearspace_consts_.xmcof = 0.0;
        if (Eccentricity() > 1.0e-4)
            nearspace_consts_.xmcof = -kTWOTHIRD * coef * BStar() * kAE / eeta;

        nearspace_consts_.delmo = pow(1.0 + common_consts_.eta * (cos(MeanAnomoly())), 3.0);
        nearspace_consts_.sinmo = sin(MeanAnomoly());

        if (!use_simple_model_) {
            const double c1sq = common_consts_.c1 * common_consts_.c1;
            nearspace_consts_.d2 = 4.0 * RecoveredSemiMajorAxis() * tsi * c1sq;
            const double temp = nearspace_consts_.d2 * tsi * common_consts_.c1 / 3.0;
            nearspace_consts_.d3 = (17.0 * RecoveredSemiMajorAxis() + s4) * temp;
            nearspace_consts_.d4 = 0.5 * temp * RecoveredSemiMajorAxis() *
                    tsi * (221.0 * RecoveredSemiMajorAxis() + 31.0 * s4) * common_consts_.c1;
            nearspace_consts_.t3cof = nearspace_consts_.d2 + 2.0 * c1sq;
            nearspace_consts_.t4cof = 0.25 * (3.0 * nearspace_consts_.d3 + common_consts_.c1 *
                    (12.0 * nearspace_consts_.d2 + 10.0 * c1sq));
            nearspace_consts_.t5cof = 0.2 * (3.0 * nearspace_consts_.d4 + 12.0 * common_consts_.c1 *
                    nearspace_consts_.d3 + 6.0 * nearspace_consts_.d2 * nearspace_consts_.d2 + 15.0 *
                    c1sq * (2.0 * nearspace_consts_.d2 + c1sq));
        }
    }

    first_run_ = false;
}

void SGP4::FindPosition(Eci* eci, double tsince) const {

    if (use_deep_space_)
        FindPositionSDP4(eci, tsince);
    else
        FindPositionSGP4(eci, tsince);
}

void SGP4::FindPosition(Eci* eci, const Julian& date) const {

    Timespan diff = date - Epoch();

    FindPosition(eci, diff.GetTotalMinutes());
}

void SGP4::FindPositionSDP4(Eci* eci, double tsince) const {

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
    double xmdf = MeanAnomoly() + common_consts_.xmdot * tsince;
    double omgadf = ArgumentPerigee() + common_consts_.omgdot * tsince;
    const double xnoddf = AscendingNode() + common_consts_.xnodot * tsince;

    const double tsq = tsince * tsince;
    xnode = xnoddf + common_consts_.xnodcf * tsq;
    double tempa = 1.0 - common_consts_.c1 * tsince;
    double tempe = BStar() * common_consts_.c4 * tsince;
    double templ = common_consts_.t2cof * tsq;

    double xn = RecoveredMeanMotion();
    e = Eccentricity();
    xincl = Inclination();

    DeepSpaceSecular(tsince, &xmdf, &omgadf, &xnode, &e, &xincl, &xn);

    if (xn <= 0.0) {
        throw SatelliteException("Error: #2 (xn <= 0.0)");
    }

    a = pow(kXKE / xn, kTWOTHIRD) * pow(tempa, 2.0);
    e = e - tempe;

    /*
     * fix tolerance for error recognition
     */
    if (e >= 1.0 || e < -0.001) {
        throw SatelliteException("Error: #1 (e >= 1.0 || e < -0.001)");
    }
    /*
     * fix tolerance to avoid a divide by zero
     */
    if (e < 1.0e-6)
        e = 1.0e-6;

    /*
    xmdf += RecoveredMeanMotion() * templ;
    double xlm = xmdf + omgadf + xnode;
    xnode = fmod(xnode, TWOPI);
    omgadf = fmod(omgadf, TWOPI);
    xlm = fmod(xlm, TWOPI);
    double xmam = fmod(xlm - omgadf - xnode, TWOPI);
     */

    double xmam = xmdf + RecoveredMeanMotion() * templ;

    DeepSpacePeriodics(tsince, &e, &xincl, &omgadf, &xnode, &xmam);

    /*
     * keeping xincl positive important unless you need to display xincl
     * and dislike negative inclinations
     */
    if (xincl < 0.0) {
        xincl = -xincl;
        xnode += kPI;
        omgadf = omgadf - kPI;
    }

    xl = xmam + omgadf + xnode;
    omega = omgadf;

    if (e < 0.0 || e > 1.0) {
        throw SatelliteException("Error: #3 (e < 0.0 || e > 1.0)");
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
        perturbed_xlcof = 0.125 * common_consts_.a3ovk2 * perturbed_sinio * (3.0 + 5.0 * perturbed_cosio) / (1.0 + perturbed_cosio);
    else
        perturbed_xlcof = 0.125 * common_consts_.a3ovk2 * perturbed_sinio * (3.0 + 5.0 * perturbed_cosio) / 1.5e-12;

    const double perturbed_aycof = 0.25 * common_consts_.a3ovk2 * perturbed_sinio;

    /*
     * using calculated values, find position and velocity
     */
    CalculateFinalPositionVelocity(eci, tsince, e,
            a, omega, xl, xnode,
            xincl, perturbed_xlcof, perturbed_aycof,
            perturbed_x3thm1, perturbed_x1mth2, perturbed_x7thm1,
            perturbed_cosio, perturbed_sinio);

}

void SGP4::FindPositionSGP4(Eci* eci, double tsince) const {

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
    const double xmdf = MeanAnomoly() + common_consts_.xmdot * tsince;
    const double omgadf = ArgumentPerigee() + common_consts_.omgdot * tsince;
    const double xnoddf = AscendingNode() + common_consts_.xnodot * tsince;

    const double tsq = tsince * tsince;
    xnode = xnoddf + common_consts_.xnodcf * tsq;
    double tempa = 1.0 - common_consts_.c1 * tsince;
    double tempe = BStar() * common_consts_.c4 * tsince;
    double templ = common_consts_.t2cof * tsq;

    xincl = Inclination();
    omega = omgadf;
    double xmp = xmdf;

    if (!use_simple_model_) {
        const double delomg = nearspace_consts_.omgcof * tsince;
        const double delm = nearspace_consts_.xmcof * (pow(1.0 + common_consts_.eta * cos(xmdf), 3.0) - nearspace_consts_.delmo);
        const double temp = delomg + delm;

        xmp += temp;
        omega = omega - temp;

        const double tcube = tsq * tsince;
        const double tfour = tsince * tcube;

        tempa = tempa - nearspace_consts_.d2 * tsq - nearspace_consts_.d3 * tcube - nearspace_consts_.d4 * tfour;
        tempe += BStar() * nearspace_consts_.c5 * (sin(xmp) - nearspace_consts_.sinmo);
        templ += nearspace_consts_.t3cof * tcube + tfour * (nearspace_consts_.t4cof + tsince * nearspace_consts_.t5cof);
    }

    a = RecoveredSemiMajorAxis() * pow(tempa, 2.0);
    e = Eccentricity() - tempe;
    xl = xmp + omega + xnode + RecoveredMeanMotion() * templ;

    if (xl <= 0.0) {
        throw SatelliteException("Error: #2 (xl <= 0.0)");
    }

    /*
     * fix tolerance for error recognition
     */
    if (e >= 1.0 || e < -0.001) {
        throw SatelliteException("Error: #1 (e >= 1.0 || e < -0.001)");
    }
    /*
     * fix tolerance to avoid a divide by zero
     */
    if (e < 1.0e-6)
        e = 1.0e-6;

    /*
     * using calculated values, find position and velocity
     * we can pass in constants from Initialize() as these dont change
     */
    CalculateFinalPositionVelocity(eci, tsince, e,
            a, omega, xl, xnode,
            xincl, common_consts_.xlcof, common_consts_.aycof,
            common_consts_.x3thm1, common_consts_.x1mth2, common_consts_.x7thm1,
            common_consts_.cosio, common_consts_.sinio);

}

void SGP4::CalculateFinalPositionVelocity(Eci* eci, const double& tsince, const double& e,
        const double& a, const double& omega, const double& xl, const double& xnode,
        const double& xincl, const double& xlcof, const double& aycof,
        const double& x3thm1, const double& x1mth2, const double& x7thm1,
        const double& cosio, const double& sinio) const {

    double temp;
    double temp1;
    double temp2;
    double temp3;

    const double beta = sqrt(1.0 - e * e);
    const double xn = kXKE / pow(a, 1.5);
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

    if (pl < 0.0) {
        throw SatelliteException("Error: #4 (pl < 0.0)");
    }

    const double r = a * (1.0 - ecose);
    temp1 = 1.0 / r;
    const double rdot = kXKE * sqrt(a) * esine * temp1;
    const double rfdot = kXKE * sqrt(pl) * temp1;
    temp2 = a * temp1;
    const double betal = sqrt(temp);
    temp3 = 1.0 / (1.0 + betal);
    const double cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
    const double sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
    const double u = atan2(sinu, cosu);
    const double sin2u = 2.0 * sinu * cosu;
    const double cos2u = 2.0 * cosu * cosu - 1.0;
    temp = 1.0 / pl;
    temp1 = kCK2 * temp;
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

    if (rk < 1.0) {
        throw SatelliteException("Error: #6 Satellite decayed (rk < 1.0)");
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
    const double x = rk * ux * kXKMPER;
    const double y = rk * uy * kXKMPER;
    const double z = rk * uz * kXKMPER;
    Vector position(x, y, z);
    const double xdot = (rdotk * ux + rfdotk * vx) * kXKMPER / 60.0;
    const double ydot = (rdotk * uy + rfdotk * vy) * kXKMPER / 60.0;
    const double zdot = (rdotk * uz + rfdotk * vz) * kXKMPER / 60.0;
    Vector velocity(xdot, ydot, zdot);

    Timespan diff;
    diff.AddMinutes(tsince);
    Julian julian = Epoch() + diff;
    (*eci) = Eci(julian, position, velocity);
}

/*
 * deep space initialization
 */
void SGP4::DeepSpaceInitialize(const double& eosq, const double& sinio, const double& cosio, const double& betao,
        const double& theta2, const double& betao2,
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

    const double aqnv = 1.0 / RecoveredSemiMajorAxis();
    const double xpidot = omgdot + xnodot;
    const double sinq = sin(AscendingNode());
    const double cosq = cos(AscendingNode());
    const double sing = sin(ArgumentPerigee());
    const double cosg = cos(ArgumentPerigee());

    /*
     * initialize lunar / solar terms
     */
    const double jday = Epoch().FromJan1_12h_1900();

    const double xnodce = 4.5236020 - 9.2422029e-4 * jday;
    const double xnodce_temp = fmod(xnodce, kTWOPI);
    const double stem = sin(xnodce_temp);
    const double ctem = cos(xnodce_temp);
    const double zcosil = 0.91375164 - 0.03568096 * ctem;
    const double zsinil = sqrt(1.0 - zcosil * zcosil);
    const double zsinhl = 0.089683511 * stem / zsinil;
    const double zcoshl = sqrt(1.0 - zsinhl * zsinhl);
    const double c = 4.7199672 + 0.22997150 * jday;
    const double gam = 5.8351514 + 0.0019443680 * jday;
    deepspace_consts_.zmol = Fmod2p(c - gam);
    double zx = 0.39785416 * stem / zsinil;
    double zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    zx = atan2(zx, zy);
    zx = fmod(gam + zx - xnodce, kTWOPI);

    const double zcosgl = cos(zx);
    const double zsingl = sin(zx);
    deepspace_consts_.zmos = Fmod2p(6.2565837 + 0.017201977 * jday);

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
        if (Inclination() < 5.2359877e-2 || Inclination() > kPI - 5.2359877e-2) {
            shdq = 0.0;
        } else {
            shdq = (-zn * s2 * (z21 + z23)) / sinio;
        }

        deepspace_consts_.ee2 = 2.0 * s1 * s6;
        deepspace_consts_.e3 = 2.0 * s1 * s7;
        deepspace_consts_.xi2 = 2.0 * s2 * z12;
        deepspace_consts_.xi3 = 2.0 * s2 * (z13 - z11);
        deepspace_consts_.xl2 = -2.0 * s3 * z2;
        deepspace_consts_.xl3 = -2.0 * s3 * (z3 - z1);
        deepspace_consts_.xl4 = -2.0 * s3 * (-21.0 - 9.0 * eosq) * ze;
        deepspace_consts_.xgh2 = 2.0 * s4 * z32;
        deepspace_consts_.xgh3 = 2.0 * s4 * (z33 - z31);
        deepspace_consts_.xgh4 = -18.0 * s4 * ze;
        deepspace_consts_.xh2 = -2.0 * s2 * z22;
        deepspace_consts_.xh3 = -2.0 * s2 * (z23 - z21);

        if (cnt == 1)
            break;
        /*
         * do lunar terms
         */
        deepspace_consts_.sse = se;
        deepspace_consts_.ssi = si;
        deepspace_consts_.ssl = sl;
        deepspace_consts_.ssh = shdq;
        deepspace_consts_.ssg = sgh - cosio * deepspace_consts_.ssh;
        deepspace_consts_.se2 = deepspace_consts_.ee2;
        deepspace_consts_.si2 = deepspace_consts_.xi2;
        deepspace_consts_.sl2 = deepspace_consts_.xl2;
        deepspace_consts_.sgh2 = deepspace_consts_.xgh2;
        deepspace_consts_.sh2 = deepspace_consts_.xh2;
        deepspace_consts_.se3 = deepspace_consts_.e3;
        deepspace_consts_.si3 = deepspace_consts_.xi3;
        deepspace_consts_.sl3 = deepspace_consts_.xl3;
        deepspace_consts_.sgh3 = deepspace_consts_.xgh3;
        deepspace_consts_.sh3 = deepspace_consts_.xh3;
        deepspace_consts_.sl4 = deepspace_consts_.xl4;
        deepspace_consts_.sgh4 = deepspace_consts_.xgh4;
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

    deepspace_consts_.sse += se;
    deepspace_consts_.ssi += si;
    deepspace_consts_.ssl += sl;
    deepspace_consts_.ssg += sgh - cosio * shdq;
    deepspace_consts_.ssh += shdq;

    deepspace_consts_.resonance_flag = false;
    deepspace_consts_.synchronous_flag = false;
    bool initialize_integrator = true;

    if (RecoveredMeanMotion() < 0.0052359877 && RecoveredMeanMotion() > 0.0034906585) {

        /*
         * 24h synchronous resonance terms initialization
         */
        deepspace_consts_.resonance_flag = true;
        deepspace_consts_.synchronous_flag = true;

        const double g200 = 1.0 + eosq * (-2.5 + 0.8125 * eosq);
        const double g310 = 1.0 + 2.0 * eosq;
        const double g300 = 1.0 + eosq * (-6.0 + 6.60937 * eosq);
        const double f220 = 0.75 * (1.0 + cosio) * (1.0 + cosio);
        const double f311 = 0.9375 * sinio * sinio * (1.0 + 3.0 * cosio) - 0.75 * (1.0 + cosio);
        double f330 = 1.0 + cosio;
        f330 = 1.875 * f330 * f330 * f330;
        deepspace_consts_.del1 = 3.0 * RecoveredMeanMotion() * RecoveredMeanMotion() * aqnv * aqnv;
        deepspace_consts_.del2 = 2.0 * deepspace_consts_.del1 * f220 * g200 * Q22;
        deepspace_consts_.del3 = 3.0 * deepspace_consts_.del1 * f330 * g300 * Q33 * aqnv;
        deepspace_consts_.del1 = deepspace_consts_.del1 * f311 * g310 * Q31 * aqnv;

        integrator_consts_.xlamo = MeanAnomoly() + AscendingNode() + ArgumentPerigee() - deepspace_consts_.gsto;
        bfact = xmdot + xpidot - kTHDT;
        bfact += deepspace_consts_.ssl + deepspace_consts_.ssg + deepspace_consts_.ssh;

    } else if (RecoveredMeanMotion() < 8.26e-3 || RecoveredMeanMotion() > 9.24e-3 || Eccentricity() < 0.5) {
        initialize_integrator = false;
    } else {
        /*
         * geopotential resonance initialization for 12 hour orbits
         */
        deepspace_consts_.resonance_flag = true;

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
            g520 = -532.114 + 3017.977 * Eccentricity() - 5740.032 * eosq + 3708.276 * eoc;
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
        deepspace_consts_.d2201 = temp * f220 * g201;
        deepspace_consts_.d2211 = temp * f221 * g211;
        temp1 = temp1 * aqnv;
        temp = temp1 * ROOT32;
        deepspace_consts_.d3210 = temp * f321 * g310;
        deepspace_consts_.d3222 = temp * f322 * g322;
        temp1 = temp1 * aqnv;
        temp = 2.0 * temp1 * ROOT44;
        deepspace_consts_.d4410 = temp * f441 * g410;
        deepspace_consts_.d4422 = temp * f442 * g422;
        temp1 = temp1 * aqnv;
        temp = temp1 * ROOT52;
        deepspace_consts_.d5220 = temp * f522 * g520;
        deepspace_consts_.d5232 = temp * f523 * g532;
        temp = 2.0 * temp1 * ROOT54;
        deepspace_consts_.d5421 = temp * f542 * g521;
        deepspace_consts_.d5433 = temp * f543 * g533;

        integrator_consts_.xlamo = MeanAnomoly() + AscendingNode() + AscendingNode() - deepspace_consts_.gsto - deepspace_consts_.gsto;
        bfact = xmdot + xnodot + xnodot - kTHDT - kTHDT;
        bfact = bfact + deepspace_consts_.ssl + deepspace_consts_.ssh + deepspace_consts_.ssh;
    }

    if (initialize_integrator) {
        /*
         * initialize integrator
         */
        integrator_consts_.xfact = bfact - RecoveredMeanMotion();
        integrator_params_.atime = 0.0;
        integrator_params_.xni = RecoveredMeanMotion();
        integrator_params_.xli = integrator_consts_.xlamo;
        /*
         * precompute dot terms for epoch
         */
        DeepSpaceCalcDotTerms(&integrator_consts_.xndot_0, &integrator_consts_.xnddt_0, &integrator_consts_.xldot_0);
    }
}

void SGP4::DeepSpaceCalculateLunarSolarTerms(const double t, double* pe, double* pinc,
        double* pl, double* pgh, double* ph) const {

    static const double ZES = 0.01675;
    static const double ZNS = 1.19459E-5;
    static const double ZNL = 1.5835218E-4;
    static const double ZEL = 0.05490;

    /*
     * calculate solar terms for time t
     */
    double zm = deepspace_consts_.zmos + ZNS * t;
    if (first_run_)
        zm = deepspace_consts_.zmos;
    double zf = zm + 2.0 * ZES * sin(zm);
    double sinzf = sin(zf);
    double f2 = 0.5 * sinzf * sinzf - 0.25;
    double f3 = -0.5 * sinzf * cos(zf);
    const double ses = deepspace_consts_.se2 * f2 + deepspace_consts_.se3 * f3;
    const double sis = deepspace_consts_.si2 * f2 + deepspace_consts_.si3 * f3;
    const double sls = deepspace_consts_.sl2 * f2 + deepspace_consts_.sl3 * f3 + deepspace_consts_.sl4 * sinzf;
    const double sghs = deepspace_consts_.sgh2 * f2 + deepspace_consts_.sgh3 * f3 + deepspace_consts_.sgh4 * sinzf;
    const double shs = deepspace_consts_.sh2 * f2 + deepspace_consts_.sh3 * f3;

    /*
     * calculate lunar terms for time t
     */
    zm = deepspace_consts_.zmol + ZNL * t;
    if (first_run_)
        zm = deepspace_consts_.zmol;
    zf = zm + 2.0 * ZEL * sin(zm);
    sinzf = sin(zf);
    f2 = 0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * cos(zf);
    const double sel = deepspace_consts_.ee2 * f2 + deepspace_consts_.e3 * f3;
    const double sil = deepspace_consts_.xi2 * f2 + deepspace_consts_.xi3 * f3;
    const double sll = deepspace_consts_.xl2 * f2 + deepspace_consts_.xl3 * f3 + deepspace_consts_.xl4 * sinzf;
    const double sghl = deepspace_consts_.xgh2 * f2 + deepspace_consts_.xgh3 * f3 + deepspace_consts_.xgh4 * sinzf;
    const double shl = deepspace_consts_.xh2 * f2 + deepspace_consts_.xh3 * f3;

    /*
     * merge calculated values
     */
    (*pe) = ses + sel;
    (*pinc) = sis + sil;
    (*pl) = sls + sll;
    (*pgh) = sghs + sghl;
    (*ph) = shs + shl;
}

/*
 * calculate lunar / solar periodics and apply
 */
void SGP4::DeepSpacePeriodics(const double& t, double* em,
        double* xinc, double* omgasm, double* xnodes, double* xll) const {

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
    DeepSpaceCalculateLunarSolarTerms(t, &pe, &pinc, &pl, &pgh, &ph);

    if (!first_run_) {

        (*xinc) += pinc;
        (*em) += pe;

        /* Spacetrack report #3 has sin/cos from before perturbations
         * added to xinc (oldxinc), but apparently report # 6 has then
         * from after they are added.
         * use for strn3
         * if (Inclination() >= 0.2)
         * use for gsfc
         * if (xinc >= 0.2)
         * (moved from start of function)
         */
        const double sinis = sin((*xinc));
        const double cosis = cos((*xinc));

        if ((*xinc) >= 0.2) {
            /*
             * apply periodics directly
             */
            const double tmp_ph = ph / sinis;

            (*omgasm) += pgh - cosis * tmp_ph;
            (*xnodes) += tmp_ph;
            (*xll) += pl;
        } else {
            /*
             * apply periodics with lyddane modification
             */
            const double sinok = sin((*xnodes));
            const double cosok = cos((*xnodes));
            double alfdp = sinis * sinok;
            double betdp = sinis * cosok;
            const double dalf = ph * cosok + pinc * cosis * sinok;
            const double dbet = -ph * sinok + pinc * cosis * cosok;

            alfdp += dalf;
            betdp += dbet;

            (*xnodes) = fmod((*xnodes), kTWOPI);
            if ((*xnodes) < 0.0)
                (*xnodes) += kTWOPI;

            double xls = (*xll) + (*omgasm) + cosis * (*xnodes);
            double dls = pl + pgh - pinc * (*xnodes) * sinis;
            xls += dls;

            /*
             * save old xnodes value
             */
            const double oldxnodes = (*xnodes);

            (*xnodes) = atan2(alfdp, betdp);
            if ((*xnodes) < 0.0)
                (*xnodes) += kTWOPI;

            /*
             * Get perturbed xnodes in to same quadrant as original.
             * RAAN is in the range of 0 to 360 degrees
             * atan2 is in the range of -180 to 180 degrees
             */
            if (fabs(oldxnodes - (*xnodes)) > kPI) {
                if ((*xnodes) < oldxnodes)
                    (*xnodes) += kTWOPI;
                else
                    (*xnodes) = (*xnodes) - kTWOPI;
            }

            (*xll) += pl;
            (*omgasm) = xls - (*xll) - cosis * (*xnodes);
        }
    }
}

/*
 * deep space secular effects
 */
void SGP4::DeepSpaceSecular(const double& t, double* xll, double* omgasm,
        double* xnodes, double* em, double* xinc, double* xn) const {

    static const double STEP = 720.0;
    static const double STEP2 = 259200.0;

    (*xll) += deepspace_consts_.ssl * t;
    (*omgasm) += deepspace_consts_.ssg * t;
    (*xnodes) += deepspace_consts_.ssh * t;
    (*em) += deepspace_consts_.sse * t;
    (*xinc) += deepspace_consts_.ssi * t;

    if (!deepspace_consts_.resonance_flag)
        return;

    /*
     * 1st condition (if t is less than one time step from epoch)
     * 2nd condition (if integrator_params_.atime and t are of opposite signs, so zero crossing required)
     * 3rd condition (if t is closer to zero than integrator_params_.atime, only integrate away from zero)
     */
    if (fabs(t) < STEP ||
            t * integrator_params_.atime <= 0.0 ||
            fabs(t) < fabs(integrator_params_.atime)) {
        /*
         * restart from epoch
         */
        integrator_params_.atime = 0.0;
        integrator_params_.xni = RecoveredMeanMotion();
        integrator_params_.xli = integrator_consts_.xlamo;

        /*
         * restore precomputed values for epoch
         */
        integrator_params_.xndot_t = integrator_consts_.xndot_0;
        integrator_params_.xnddt_t = integrator_consts_.xnddt_0;
        integrator_params_.xldot_t = integrator_consts_.xldot_0;
    }

    double ft = t - integrator_params_.atime;

    /*
     * if time difference (ft) is greater than the time step (720.0)
     * loop around until integrator_params_.atime is within one time step of t
     */
    if (fabs(ft) >= STEP) {

        /*
         * calculate step direction to allow integrator_params_.atime
         * to catch up with t
         */
        double delt = -STEP;
        if (ft >= 0.0)
            delt = STEP;

        do {
            /*
             * integrate using current dot terms
             */
            DeepSpaceIntegrator(delt, STEP2, integrator_params_.xndot_t, integrator_params_.xnddt_t, integrator_params_.xldot_t);

            /*
             * calculate dot terms for next integration
             */
            DeepSpaceCalcDotTerms(&integrator_params_.xndot_t, &integrator_params_.xnddt_t, &integrator_params_.xldot_t);

            ft = t - integrator_params_.atime;
        } while (fabs(ft) >= STEP);
    }

    /*
     * integrator
     */
    (*xn) = integrator_params_.xni + integrator_params_.xndot_t * ft + integrator_params_.xnddt_t * ft * ft * 0.5;
    const double xl = integrator_params_.xli + integrator_params_.xldot_t * ft + integrator_params_.xndot_t * ft * ft * 0.5;
    const double temp = -(*xnodes) + deepspace_consts_.gsto + t * kTHDT;

    if (deepspace_consts_.synchronous_flag)
        (*xll) = xl + temp - (*omgasm);
    else
        (*xll) = xl + temp + temp;
}

/*
 * calculate dot terms
 */
void SGP4::DeepSpaceCalcDotTerms(double* xndot, double* xnddt, double* xldot) const {

    static const double G22 = 5.7686396;
    static const double G32 = 0.95240898;
    static const double G44 = 1.8014998;
    static const double G52 = 1.0508330;
    static const double G54 = 4.4108898;
    static const double FASX2 = 0.13130908;
    static const double FASX4 = 2.8843198;
    static const double FASX6 = 0.37448087;

    if (deepspace_consts_.synchronous_flag) {

        (*xndot) = deepspace_consts_.del1 * sin(integrator_params_.xli - FASX2) +
                deepspace_consts_.del2 * sin(2.0 * (integrator_params_.xli - FASX4)) +
                deepspace_consts_.del3 * sin(3.0 * (integrator_params_.xli - FASX6));
        (*xnddt) = deepspace_consts_.del1 * cos(integrator_params_.xli - FASX2) + 2.0 *
                deepspace_consts_.del2 * cos(2.0 * (integrator_params_.xli - FASX4)) + 3.0 *
                deepspace_consts_.del3 * cos(3.0 * (integrator_params_.xli - FASX6));

    } else {

        const double xomi = ArgumentPerigee() + common_consts_.omgdot * integrator_params_.atime;
        const double x2omi = xomi + xomi;
        const double x2li = integrator_params_.xli + integrator_params_.xli;

        (*xndot) = deepspace_consts_.d2201 * sin(x2omi + integrator_params_.xli - G22)
                + deepspace_consts_.d2211 * sin(integrator_params_.xli - G22)
                + deepspace_consts_.d3210 * sin(xomi + integrator_params_.xli - G32)
                + deepspace_consts_.d3222 * sin(-xomi + integrator_params_.xli - G32)
                + deepspace_consts_.d4410 * sin(x2omi + x2li - G44)
                + deepspace_consts_.d4422 * sin(x2li - G44)
                + deepspace_consts_.d5220 * sin(xomi + integrator_params_.xli - G52)
                + deepspace_consts_.d5232 * sin(-xomi + integrator_params_.xli - G52)
                + deepspace_consts_.d5421 * sin(xomi + x2li - G54)
                + deepspace_consts_.d5433 * sin(-xomi + x2li - G54);
        (*xnddt) = deepspace_consts_.d2201 * cos(x2omi + integrator_params_.xli - G22)
                + deepspace_consts_.d2211 * cos(integrator_params_.xli - G22)
                + deepspace_consts_.d3210 * cos(xomi + integrator_params_.xli - G32)
                + deepspace_consts_.d3222 * cos(-xomi + integrator_params_.xli - G32)
                + deepspace_consts_.d5220 * cos(xomi + integrator_params_.xli - G52)
                + deepspace_consts_.d5232 * cos(-xomi + integrator_params_.xli - G52)
                + 2.0 * (deepspace_consts_.d4410 * cos(x2omi + x2li - G44)
                + deepspace_consts_.d4422 * cos(x2li - G44)
                + deepspace_consts_.d5421 * cos(xomi + x2li - G54)
                + deepspace_consts_.d5433 * cos(-xomi + x2li - G54));
    }

    (*xldot) = integrator_params_.xni + integrator_consts_.xfact;
    (*xnddt) = (*xnddt) * (*xldot);
}

/*
 * deep space integrator for time period of delt
 */
void SGP4::DeepSpaceIntegrator(const double delt, const double step2,
        const double xndot, const double xnddt, const double xldot) const {

    /*
     * integrator
     */
    integrator_params_.xli += xldot * delt + xndot * step2;
    integrator_params_.xni += xndot * delt + xnddt * step2;

    /*
     * increment integrator time
     */
    integrator_params_.atime += delt;
}

void SGP4::Reset() {

    /*
     * common variables
     */
    first_run_ = true;
    use_simple_model_ = false;
    use_deep_space_ = false;

    memset(&common_consts_, 0, sizeof (common_consts_));
    memset(&nearspace_consts_, 0, sizeof (nearspace_consts_));
    memset(&deepspace_consts_, 0, sizeof (deepspace_consts_));
    memset(&integrator_consts_, 0, sizeof (integrator_consts_));
    memset(&integrator_params_, 0, sizeof (integrator_params_));

    mean_anomoly_ = ascending_node_ = argument_perigee_ = eccentricity_ =
            inclination_ = mean_motion_ = bstar_ = recovered_semi_major_axis_ =
            recovered_mean_motion_ = perigee_ = period_ = 0.0;

    epoch_ = Julian();
}
