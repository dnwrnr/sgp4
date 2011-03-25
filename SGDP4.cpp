#include "SGDP4.h"

#include "SatelliteException.h"

#include <math.h>

SGDP4::SGDP4(void) {
    first_run_ = true;
}

SGDP4::~SGDP4(void) {
}

void SGDP4::Initialize(const Tle& tle) {
    /*
     * extract and format tle data
     */
    tle_data_.bstar = tle.GetField(Tle::FLD_BSTAR);
    tle_data_.eo = tle.GetField(Tle::FLD_E);
    tle_data_.omega = tle.GetField(Tle::FLD_ARGPER, Tle::U_RAD);
    tle_data_.xincl = tle.GetField(Tle::FLD_I, Tle::U_RAD);
    tle_data_.xmo = tle.GetField(Tle::FLD_M, Tle::U_RAD);
    tle_data_.xno = tle.GetField(Tle::FLD_MMOTION) / (1440.0 / Globals::TWOPI());
    tle_data_.xnodeo = tle.GetField(Tle::FLD_RAAN, Tle::U_RAD);

    /*
     * generate julian date for tle epoch
     */
    int year = static_cast<int> (tle.GetField(Tle::FLD_EPOCHYEAR));
    if (year < 57)
        year += 2000;
    else
        year += 1900;
    double day = tle.GetField(Tle::FLD_EPOCHDAY);
    Julian jul(year, day);
    tle_data_.epoch = jul;

    double cosio_ = 0.0;
    double sinio_ = 0.0;
    double betao2_ = 0.0;
    double betao_ = 0.0;
    double xnodp_ = 0.0;
    double aodp_ = 0.0;

    double gsto = 0.0;

    bool use_simple_model_ = false;
    bool use_deep_space_ = false;

    /*
     * recover original mean motion (xnodp) and semimajor axis (aodp)
     * from input elements
     */
    double a1 = pow(Globals::XKE() / tle_data_.xno, Globals::TOTHRD());
    cosio_ = cos(tle_data_.xincl);
    sinio_ = sin(tle_data_.xincl);
    double theta2 = cosio_ * cosio_;
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = tle_data_.eo * tle_data_.eo;
    betao2_ = 1.0 - eosq;
    betao_ = sqrt(betao2_);
    double del1 = 1.5 * Globals::CK2() * x3thm1 / (a1 * a1 * betao_ * betao2_);
    double ao = a1 * (1.0 - del1 * (0.5 * Globals::TOTHRD() + del1 * (1.0 + 134.0 / 81.0 * del1)));
    double delo = 1.5 * Globals::CK2() * x3thm1 / (ao * ao * betao_ * betao2_);
    /*
     * recovered mean motion
     */
    xnodp_ = tle_data_.xno / (1.0 + delo);
    /*
     * recovered semimajor axis
     */
    aodp_ = ao / (1.0 - delo);

    gsto = tle_data_.epoch.ToGMST();

    double rp = aodp_ * (1.0 - tle_data_.eo);
    double perigee = (aodp_ * (1.0 - tle_data_.eo) - Globals::AE()) * Globals::XKMPER();
    double period = Globals::TWOPI() / xnodp_;

    /*
     * for perigee less than 220 kilometers, the simple_model flag is set and
     * the equations are truncated to linear variation in sqrt a and
     * quadratic variation in mean anomly. also, the c3 term, the
     * delta omega term and the delta m term are dropped
     */
    use_simple_model_ = false;
    if (rp < (220.0 / Globals::XKMPER() + Globals::AE()))
        use_simple_model_ = true;

    double s4_ = 0.0;
    double qoms24_ = 0.0;

    s4_ = Globals::S();
    qoms24_ = Globals::QOMS2T();
    /*
     * for perigee below 156km, the values of
     * s4 and qoms2t are altered
     */
    if (perigee < 156.0) {
        s4_ = perigee - 78.0;
        if (perigee <= 98.0) {
            s4_ = 20.0;
        }

        qoms24_ = pow((120.0 - s4_) * Globals::AE() / Globals::XKMPER(), 4.0);
        s4_ = s4_ / Globals::XKMPER() + Globals::AE();
    }

    if (period >= 225.0) {
        use_deep_space_ = true;
    }

    double pinvsq = 1.0 / (aodp_ * aodp_ * betao2_ * betao2_);
    double sing = sin(tle_data_.omega);
    double cosg = cos(tle_data_.omega);

    double tsi_ = 0.0;
    double eta_ = 0.0;
    double eeta_ = 0.0;
    double coef_ = 0.0;
    double coef1_ = 0.0;
    double c1_ = 0.0;
    double a3ovk2_ = 0.0;

    tsi_ = 1.0 / (aodp_ - s4_);
    eta_ = aodp_ * tle_data_.eo * tsi_;
    double etasq = eta_ * eta_;
    eeta_ = tle_data_.eo * eta_;
    double psisq = fabs(1.0 - etasq);
    coef_ = qoms24_ * pow(tsi_, 4.0);
    coef1_ = coef_ / pow(psisq, 3.5);
    double c2 = coef1_ * xnodp_ * (aodp_ * (1.0 + 1.5 * etasq + eeta_ *
            (4.0 + etasq)) + 0.75 * Globals::CK2() * tsi_ / psisq *
            x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    c1_ = tle_data_.bstar * c2;
    a3ovk2_ = -Globals::XJ3() / Globals::CK2() * pow(Globals::AE(), 3.0);

    double c3_ = 0.0;
    double c4_ = 0.0;
    double c5_ = 0.0;
    double xmdot_ = 0.0;
    double omgdot_ = 0.0;
    double xnodot_ = 0.0;
    double xnodcf_ = 0.0;
    double t2cof_ = 0.0;
    double xlcof_ = 0.0; // move to end?
    double aycof_ = 0.0; // move to end?
    double x7thm1_ = 0.0; // move to end?
    double omgcof_ = 0.0;
    double xmcof_ = 0.0;
    double delmo_ = 0.0;
    double sinmo_ = 0.0;

    double x1mth2 = 1.0 - theta2;
    c4_ = 2.0 * xnodp_ * coef1_ * aodp_ * betao2_ *
            (eta_ * (2.0 + 0.5 * etasq) + tle_data_.eo * (0.5 + 2.0 * etasq) -
            2.0 * Globals::CK2() * tsi_ / (aodp_ * psisq) *
            (-3.0 * x3thm1 * (1.0 - 2.0 * eeta_ + etasq *
            (1.5 - 0.5 * eeta_)) + 0.75 * x1mth2 * (2.0 * etasq - eeta_ *
            (1.0 + etasq)) * cos(2.0 * tle_data_.omega)));
    double theta4 = theta2 * theta2;
    double temp1 = 3.0 * Globals::CK2() * pinvsq * xnodp_;
    double temp2 = temp1 * Globals::CK2() * pinvsq;
    double temp3 = 1.25 * Globals::CK4() * pinvsq * pinvsq * xnodp_;
    xmdot_ = xnodp_ + 0.5 * temp1 * betao_ * x3thm1 + 0.0625 * temp2 * betao_ *
            (13.0 - 78.0 * theta2 + 137.0 * theta4);
    double x1m5th = 1.0 - 5.0 * theta2;
    omgdot_ = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7.0 - 114.0 * theta2 +
            395.0 * theta4) + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
    double xhdot1_ = -temp1 * cosio_;
    xnodot_ = xhdot1_ + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 *
            (3.0 - 7.0 * theta2)) * cosio_;
    xnodcf_ = 3.5 * betao2_ * xhdot1_ * c1_;
    t2cof_ = 1.5 * c1_;

    if (fabs(cosio_ + 1.0) > 1.5e-12)
        xlcof_ = 0.125 * a3ovk2_ * sinio_ * (3.0 + 5.0 * cosio_) / (1.0 + cosio_);
    else
        xlcof_ = 0.125 * a3ovk2_ * sinio_ * (3.0 + 5.0 * cosio_) / 1.5e-12;

    aycof_ = 0.25 * a3ovk2_ * sinio_;
    x7thm1_ = 7.0 * theta2 - 1.0;

    if (!use_deep_space_) {

        // check result (different to telsko)
        c3_ = 0.0;
        if (tle_data_.eo > 1.0e-4)
            c3_ = coef_ * tsi_ * a3ovk2_ * xnodp_ * Globals::AE() *
            sinio_ / tle_data_.eo;

        c5_ = 2.0 * coef1_ * aodp_ * betao2_ * (1.0 + 2.75 * (etasq + eeta_) + eeta_ * etasq);
        omgcof_ = tle_data_.bstar * c3_ * cos(tle_data_.omega);

        xmcof_ = 0.0;
        if (tle_data_.eo > 1.0e-4)
            xmcof_ = -Globals::TOTHRD() * coef_ * tle_data_.bstar * Globals::AE() / eeta_;

        delmo_ = pow(1.0 + eta_ * (cos(tle_data_.xmo)), 3.0);
        sinmo_ = sin(tle_data_.xmo);
    }

    double d2_ = 0.0;
    double d3_ = 0.0;
    double d4_ = 0.0;
    double t3cof_ = 0.0;
    double t4cof_ = 0.0;
    double t5cof_ = 0.0;

    if (!use_simple_model_) {
        double c1sq = c1_ * c1_;
        d2_ = 4.0 * aodp_ * tsi_ * c1sq;
        double temp = d2_ * tsi_ * c1_ / 3.0;
        d3_ = (17.0 * aodp_ + s4_) * temp;
        d4_ = 0.5 * temp * aodp_ * tsi_ * (221.0 * aodp_ + 31.0 * s4_) * c1_;
        t3cof_ = d2_ + 2.0 * c1sq;
        t4cof_ = 0.25 * (3.0 * d3_ + c1_ * (12.0 * d2_ + 10.0 * c1sq));
        t5cof_ = 0.2 * (3.0 * d4_ + 12.0 * c1_ * d3_ + 6.0 * d2_ * d2_ + 15.0 *
                c1sq * (2.0 * d2_ + c1sq));
    } else if (use_deep_space_) {

    }

    first_run_ = false;
}
