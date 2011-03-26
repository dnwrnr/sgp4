#include "SGDP4.h"

#include "SatelliteException.h"

#include <math.h>

SGDP4::SGDP4(void) {
    first_run_ = true;
}

SGDP4::~SGDP4(void) {
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
    mean_motion_ = tle.GetField(Tle::FLD_MMOTION) / (1440.0 / Globals::TWOPI());
    bstar_ = tle.GetField(Tle::FLD_BSTAR);

    /*
     * generate julian date for tle epoch
     */
    epoch_ = tle.GetEpoch();

    /*
     * recover original mean motion (xnodp) and semimajor axis (aodp)
     * from input elements
     */
    double a1 = pow(Globals::XKE() / MeanMotion(), Globals::TOTHRD());
    cosio_ = cos(Inclination());
    double theta2 = cosio_ * cosio_;
    x3thm1_ = 3.0 * theta2 - 1.0;
    double eosq = Eccentricity() * Eccentricity();
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);
    double del1 = 1.5 * Globals::CK2() * x3thm1_ / (a1 * a1 * betao * betao2);
    double ao = a1 * (1.0 - del1 * (0.5 * Globals::TOTHRD() + del1 * (1.0 + 134.0 / 81.0 * del1)));
    double delo = 1.5 * Globals::CK2() * x3thm1_ / (ao * ao * betao * betao2);

    recovered_mean_motion_ = MeanMotion() / (1.0 + delo);
    recovered_semi_major_axis_ = ao / (1.0 - delo);

    /*
     * find perigee and period
     */
    perigee_ = (RecoveredSemiMajorAxis() * (1.0 - Eccentricity()) - Globals::AE()) * Globals::XKMPER();
    period_ = Globals::TWOPI() / RecoveredMeanMotion();

    Initialize(theta2, betao2, betao);
}

void SGDP4::Initialize(const double& theta2, const double& betao2, const double& betao) {

    cosio_ = 0.0;
    sinio_ = 0.0;

    eta_ = 0.0;
    coef1_ = 0.0;
    c1_ = 0.0;
    a3ovk2_ = 0.0;
    x1mth2_ = 0.0;
    c4_ = 0.0;
    c5_ = 0.0;

    xmdot_ = 0.0;
    omgdot_ = 0.0;
    xnodot_ = 0.0;
    xnodcf_ = 0.0;
    t2cof_ = 0.0;
    xlcof_ = 0.0;
    aycof_ = 0.0;
    x7thm1_ = 0.0;
    omgcof_ = 0.0;
    xmcof_ = 0.0;
    delmo_ = 0.0;
    sinmo_ = 0.0;

    d2_ = 0.0;
    d3_ = 0.0;
    d4_ = 0.0;
    t3cof_ = 0.0;
    t4cof_ = 0.0;
    t5cof_ = 0.0;

    gsto_ = 0.0;

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

    double s4_ = Globals::S();
    double qoms24_ = Globals::QOMS2T();
    /*
     * for perigee below 156km, the values of
     * s4 and qoms2t are altered
     */
    if (Perigee() < 156.0) {
        s4_ = Perigee() - 78.0;
        if (Perigee() <= 98.0) {
            s4_ = 20.0;
        }
        qoms24_ = pow((120.0 - s4_) * Globals::AE() / Globals::XKMPER(), 4.0);
        s4_ = s4_ / Globals::XKMPER() + Globals::AE();
    }

    double pinvsq = 1.0 / (RecoveredSemiMajorAxis() * RecoveredSemiMajorAxis() * betao2 * betao2);

    double tsi = 1.0 / (RecoveredSemiMajorAxis() - s4_);
    eta_ = RecoveredSemiMajorAxis() * Eccentricity() * tsi;
    double etasq = eta_ * eta_;
    double eeta = Eccentricity() * eta_;
    double psisq = fabs(1.0 - etasq);
    double coef = qoms24_ * pow(tsi, 4.0);
    coef1_ = coef / pow(psisq, 3.5);
    double c2 = coef1_ * RecoveredMeanMotion() * (RecoveredSemiMajorAxis() * (1.0 + 1.5 * etasq + eeta *
            (4.0 + etasq)) + 0.75 * Globals::CK2() * tsi / psisq *
            x3thm1_ * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    c1_ = BStar() * c2;
    a3ovk2_ = -Globals::XJ3() / Globals::CK2() * pow(Globals::AE(), 3.0);

    x1mth2_ = 1.0 - theta2;
    c4_ = 2.0 * RecoveredMeanMotion() * coef1_ * RecoveredSemiMajorAxis() * betao2 *
            (eta_ * (2.0 + 0.5 * etasq) + Eccentricity() * (0.5 + 2.0 * etasq) -
            2.0 * Globals::CK2() * tsi / (RecoveredSemiMajorAxis() * psisq) *
            (-3.0 * x3thm1_ * (1.0 - 2.0 * eeta + etasq *
            (1.5 - 0.5 * eeta)) + 0.75 * x1mth2_ * (2.0 * etasq - eeta *
            (1.0 + etasq)) * cos(2.0 * ArgumentPerigee())));
    double theta4 = theta2 * theta2;
    double temp1 = 3.0 * Globals::CK2() * pinvsq * RecoveredMeanMotion();
    double temp2 = temp1 * Globals::CK2() * pinvsq;
    double temp3 = 1.25 * Globals::CK4() * pinvsq * pinvsq * RecoveredMeanMotion();
    xmdot_ = RecoveredMeanMotion() + 0.5 * temp1 * betao * x3thm1_ + 0.0625 * temp2 * betao *
            (13.0 - 78.0 * theta2 + 137.0 * theta4);
    double x1m5th = 1.0 - 5.0 * theta2;
    omgdot_ = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7.0 - 114.0 * theta2 +
            395.0 * theta4) + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
    double xhdot1_ = -temp1 * cosio_;
    xnodot_ = xhdot1_ + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 *
            (3.0 - 7.0 * theta2)) * cosio_;
    xnodcf_ = 3.5 * betao2 * xhdot1_ * c1_;
    t2cof_ = 1.5 * c1_;

    if (fabs(cosio_ + 1.0) > 1.5e-12)
        xlcof_ = 0.125 * a3ovk2_ * sinio_ * (3.0 + 5.0 * cosio_) / (1.0 + cosio_);
    else
        xlcof_ = 0.125 * a3ovk2_ * sinio_ * (3.0 + 5.0 * cosio_) / 1.5e-12;

    aycof_ = 0.25 * a3ovk2_ * sinio_;
    x7thm1_ = 7.0 * theta2 - 1.0;

    if (!use_deep_space_) {

        double c3 = 0.0;
        if (Eccentricity() > 1.0e-4) {
            c3 = coef * tsi * a3ovk2_ * RecoveredMeanMotion() * Globals::AE() *
                    sinio_ / Eccentricity();
        }

        c5_ = 2.0 * coef1_ * RecoveredSemiMajorAxis() * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
        omgcof_ = BStar() * c3 * cos(ArgumentPerigee());

        xmcof_ = 0.0;
        if (Eccentricity() > 1.0e-4)
            xmcof_ = -Globals::TOTHRD() * coef * BStar() * Globals::AE() / eeta;

        delmo_ = pow(1.0 + eta_ * (cos(MeanAnomoly())), 3.0);
        sinmo_ = sin(MeanAnomoly());
    }

    if (!use_simple_model_) {
        double c1sq = c1_ * c1_;
        d2_ = 4.0 * RecoveredSemiMajorAxis() * tsi * c1sq;
        double temp = d2_ * tsi * c1_ / 3.0;
        d3_ = (17.0 * RecoveredSemiMajorAxis() + s4_) * temp;
        d4_ = 0.5 * temp * RecoveredSemiMajorAxis() * tsi * (221.0 * RecoveredSemiMajorAxis() + 31.0 * s4_) * c1_;
        t3cof_ = d2_ + 2.0 * c1sq;
        t4cof_ = 0.25 * (3.0 * d3_ + c1_ * (12.0 * d2_ + 10.0 * c1sq));
        t5cof_ = 0.2 * (3.0 * d4_ + 12.0 * c1_ * d3_ + 6.0 * d2_ * d2_ + 15.0 *
                c1sq * (2.0 * d2_ + c1sq));
    } else if (use_deep_space_) {
        gsto_ = Epoch().ToGMST();
        double sing = sin(ArgumentPerigee());
        double cosg = cos(ArgumentPerigee());
        //DeepSpaceInitialize(eosq, sinio, cosio, betao);
        //CALL DPINIT(EOSQ,SINIO,COSIO,BETAO,AODP,THETA2,
        // SING,COSG,BETAO2,XMDOT,OMGDOT,XNODOT,XNODP)
    }

    first_run_ = false;
}

void SGDP4::FindPosition(double tsince) {

    struct TleData tle_data_tsince_;
    memset(&tle_data_tsince_, 0, sizeof (tle_data_tsince_));

    tle_data_tsince_.bstar = BStar();
    tle_data_tsince_.eo = Eccentricity();
    tle_data_tsince_.omega = ArgumentPerigee();
    tle_data_tsince_.xincl = Inclination();
    tle_data_tsince_.xmo = MeanAnomoly();
    tle_data_tsince_.xno = RecoveredMeanMotion();
    tle_data_tsince_.xnodeo = AscendingNode();
    tle_data_tsince_.epoch = Epoch();

    double xl = 0.0;
    double a = 0.0;

    /*
     * update for secular gravity and atmospheric drag
     */
    double xmdf = MeanAnomoly() + xmdot_ * tsince;
    double omgadf = ArgumentPerigee() + omgdot_ * tsince;
    double xnoddf = tle_data_tsince_.xnodeo + xnodot_ * tsince;

    double tsq = tsince * tsince;
    double xnode = xnoddf + xnodcf_ * tsq;
    double tempa = 1.0 - c1_ * tsince;
    double tempe = tle_data_tsince_.bstar * c4_ * tsince;
    double templ = t2cof_ * tsq;

    tle_data_tsince_.omega = omgadf;

    if (use_deep_space_) {
        double xn = RecoveredMeanMotion();
#if 0
        CALL DPSEC(xmdf, tle_data_tsince_.omega, XNODE, tle_data_tsince_.eo, tle_data_tsince_.xincl, xn, tsince);
#endif
        a = pow(Globals::XKE() / xn, Globals::TOTHRD()) * pow(tempa, 2.0);
        tle_data_tsince_.eo -= tempe;
        double xmam = xmdf + RecoveredMeanMotion() * templ;
#if 0
        CALL DPPER(tle_data_tsince_.eo, tle_data_tsince_.xincl, tle_data_tsince_.omega, tle_data_tsince_.xnodeo, xmam);
#endif
        xl = xmam + tle_data_tsince_.omega + xnode;

        /*
         * re-compute the perturbed values
         */
        sinio_ = sin(tle_data_tsince_.xincl);
        cosio_ = cos(tle_data_tsince_.xincl);

        double theta2 = cosio_ * cosio_;

        x3thm1_ = 3.0 * theta2 - 1.0;
        x1mth2_ = 1.0 - theta2;
        x7thm1_ = 7.0 * theta2 - 1.0;

        if (fabs(cosio_ + 1.0) > 1.5e-12)
            xlcof_ = 0.125 * a3ovk2_ * sinio_ * (3.0 + 5.0 * cosio_) / (1.0 + cosio_);
        else
            xlcof_ = 0.125 * a3ovk2_ * sinio_ * (3.0 + 5.0 * cosio_) / 1.5e-12;

        aycof_ = 0.25 * a3ovk2_ * sinio_;

    } else {
        double xmp = xmdf;
        if (!use_simple_model_) {
            double delomg = omgcof_ * tsince;
            double delm = xmcof_ * (pow(1.0 + eta_ * cos(xmdf), 3.0) - delmo_);
            double temp1 = delomg + delm;
            xmp = xmdf + temp1;
            tle_data_tsince_.omega -= temp1;
            double tcube = tsq * tsince;
            double tfour = tsince * tcube;
            tempa -= d2_ * tsq - d3_ * tcube - d4_ * tfour;
            tempe += tle_data_tsince_.bstar * c5_ * (sin(xmp) - sinmo_);
            templ += t3cof_ * tcube + tfour * (t4cof_ + tsince * t5cof_);
        }
        a = RecoveredSemiMajorAxis() * pow(tempa, 2.0);
        tle_data_tsince_.eo = Eccentricity() - tempe;
        xl = xmp + tle_data_tsince_.omega + xnode + RecoveredMeanMotion() * templ;
    }

    double beta = sqrt(1.0 - tle_data_tsince_.eo * tle_data_tsince_.eo);
    double xn = Globals::XKE() / pow(a, 1.5);
    /*
     * long period periodics
     */
    double axn;
    double xll;
    double aynl;
    double xlt;
    double ayn;
    {
        axn = tle_data_tsince_.eo * cos(tle_data_tsince_.omega);
        double temp1 = 1.0 / (a * beta * beta);
        xll = temp1 * xlcof_ * axn;
        aynl = temp1 * aycof_;
        xlt = xl + xll;
        ayn = tle_data_tsince_.eo * sin(tle_data_tsince_.omega) + aynl;
    }
    /*
     * solve keplers equation
     */
    double capu = fmod(xlt - xnode, Globals::TWOPI());
    double epw = capu;

    double sinepw = 0.0;
    double cosepw = 0.0;
    double ecose = 0.0;
    double esine = 0.0;

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
            double df = 1.0 - ecose;
            double nr = f / df;

            /*
             * Newton-Raphson correction of -F/DF
             */
            epw += nr;
        }
    }
    /*
     * short period preliminary quantities
     */
    double pl;
    double r;
    double betal;
    double rdot;
    double rfdot;
    double u;
    double sin2u;
    double cos2u;
    {
        double elsq = axn * axn + ayn * ayn;
        double temp1 = 1.0 - elsq;
        pl = a * temp1;
        r = a * (1.0 - ecose);
        betal = sqrt(temp1);
        double temp2 = 1.0 / r;
        double temp3 = a * temp2;
        double temp4 = 1.0 / (1.0 + betal);
        rdot = Globals::XKE() * sqrt(a) * esine * temp2;
        rfdot = Globals::XKE() * sqrt(pl) * temp2;
        double cosu = temp3 * (cosepw - axn + ayn * esine * temp4);
        double sinu = temp3 * (sinepw - ayn - axn * esine * temp4);
        u = atan2(sinu, cosu);
        sin2u = 2.0 * sinu * cosu;
        cos2u = 2.0 * cosu * cosu - 1.0;
    }

    /*
     * update for short periodics
     */
    double rk;
    double uk;
    double xnodek;
    double xinck;
    double rdotk;
    double rfdotk;
    {
        double temp1 = 1.0 / pl;
        double temp2 = Globals::CK2() * temp1;
        double temp3 = temp2 * temp1;
        rk = r * (1.0 - 1.5 * temp3 * betal * x3thm1_) + 0.5 * temp2 * x1mth2_ * cos2u;
        uk = u - 0.25 * temp3 * x7thm1_ * sin2u;
        xnodek = xnode + 1.5 * temp3 * cosio_ * sin2u;
        xinck = tle_data_tsince_.xincl + 1.5 * temp3 * cosio_ * sinio_ * cos2u;
        rdotk = rdot - xn * temp2 * x1mth2_ * sin2u;
        rfdotk = rfdot + xn * temp2 * (x1mth2_ * cos2u + 1.5 * x3thm1_);
    }
    /*
     * orientation vectors
     */
    double sinuk = sin(uk);
    double cosuk = cos(uk);
    double sinik = sin(xinck);
    double cosik = cos(xinck);
    double sinnok = sin(xnodek);
    double cosnok = cos(xnodek);
    double xmx = -sinnok * cosik;
    double xmy = cosnok * cosik;
    double ux = xmx * sinuk + cosnok * cosuk;
    double uy = xmy * sinuk + sinnok * cosuk;
    double uz = sinik * sinuk;
    double vx = xmx * cosuk - cosnok * sinuk;
    double vy = xmy * cosuk - sinnok * sinuk;
    double vz = sinik * cosuk;
    /*
     * position and velocity
     */
    double x = rk * ux * Globals::XKMPER();
    double y = rk * uy * Globals::XKMPER();
    double z = rk * uz * Globals::XKMPER();
    double xdot = (rdotk * ux + rfdotk * vx) * Globals::XKMPER() / 60.0;
    double ydot = (rdotk * uy + rfdotk * vy) * Globals::XKMPER() / 60.0;
    double zdot = (rdotk * uz + rfdotk * vz) * Globals::XKMPER() / 60.0;
}

/*
 * entry dpinit(const double& eqsq, const double& siniq, const double& cosiq,
 * const double& rteqsq, cosq2, sinomo, cosomo, bsq, xlldot, omgdt, xnodot, xnodp)
 */

/*
 * deep space initialization
 */
void SGDP4::DeepSpaceInitialize() {

    double omgdt = 0.0;
    double cosiq = 0.0;
    double siniq = 0.0;
    double cosomo = 0.0;
    double sinomo = 0.0;
    double eqsq = 0.0;
    double bsq = 0.0;
    double shdq = 0.0;
    double cosq2 = 0.0;
    double rteqsq = 0.0;
    double xlldot = 0.0;
    double xnodot = 0.0;




    double a1;
    double a2;
    double a3;
    double a4;
    double a5;
    double a6;
    double a7;
    double a8;
    double a9;
    double a10;

    double x1;
    double x2;
    double x3;
    double x4;
    double x5;
    double x6;
    double x7;
    double x8;

    double z1;
    double z2;
    double z3;

    double z11;
    double z12;
    double z13;

    double z21;
    double z22;
    double z23;

    double z31;
    double z32;
    double z33;

    double s1;
    double s2;
    double s3;
    double s4;
    double s5;
    double s6;
    double s7;

    double se;
    double si;
    double sl;
    double sgh;

    double bfact;

    double ZNS = 1.19459E-5;
    double C1SS = 2.9864797E-6;
    double ZES = 0.01675;
    double ZNL = 1.5835218E-4;
    double C1L = 4.7968065E-7;
    double ZEL = 0.05490;
    double ZCOSIS = 0.91744867;
    double ZSINI = 0.39785416;
    double ZSINGS = -0.98088458;
    double ZCOSGS = 0.1945905;
    double Q22 = 1.7891679E-6;
    double Q31 = 2.1460748E-6;
    double Q33 = 2.2123015E-7;
    double ROOT22 = 1.7891679E-6;
    double ROOT32 = 3.7393792E-7;
    double ROOT44 = 7.3636953E-9;
    double ROOT52 = 1.1428639E-7;
    double ROOT54 = 2.1765803E-9;
    double THDT = 4.3752691E-3;

    double aqnv = 1.0 / RecoveredSemiMajorAxis();
    double xpidot = omgdt + xnodot_;
    double sinq = sin(AscendingNode());
    double cosq = cos(AscendingNode());

    /*
     * initialize lunar / solar terms
     */
    d_day_ = Epoch().FromJan1_12h_1900();

    double xnodce = 4.5236020 - 9.2422029e-4 * d_day_;
    double stem = sin(xnodce);
    double ctem = cos(xnodce);
    double zcosil = 0.91375164 - 0.03568096 * ctem;
    double zsinil = sqrt(1.0 - zcosil * zcosil);
    double zsinhl = 0.089683511 * stem / zsinil;
    double zcoshl = sqrt(1.0 - zsinhl * zsinhl);
    double c = 4.7199672 + 0.22997150 * d_day_;
    double gam = 5.8351514 + 0.0019443680 * d_day_;
    double zmol = Globals::Fmod2p(c - gam);
    double zx = 0.39785416 * stem / zsinil;
    double zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    /*
     * todo: check
     */
    zx = atan2(zx, zy);
    zx = fmod(gam + zx - xnodce, Globals::TWOPI());

    double zcosgl = cos(zx);
    double zsingl = sin(zx);
    double zmos = 6.2565837 + 0.017201977 * d_day_;
    zmos = Globals::Fmod2p(zmos);

    /*
     * do solar terms
     */
    d_savtsn_ = 1.0e20;
    double zcosg = ZCOSGS;
    double zsing = ZSINGS;
    double zcosi = ZCOSIS;
    double zsini = ZSINI;
    double zcosh = cosq;
    double zsinh = sinq;
    double cc = C1SS;
    double zn = ZNS;
    double ze = ZES;
    double zmo = d_zmos_;
    double xnoi = 1.0 / RecoveredMeanMotion();

    for (int cnt = 0; cnt < 2; cnt++) {
        /*
         * solar terms are done a second time after lunar terms are done
         */
        a1 = zcosg * zcosh + zsing * zcosi * zsinh;
        a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
        a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
        a8 = zsing * zsini;
        a9 = zsing * zsinh + zcosg * zcosi*zcosh;
        a10 = zcosg * zsini;
        a2 = cosiq * a7 + siniq * a8;
        a4 = cosiq * a9 + siniq * a10;
        a5 = -siniq * a7 + cosiq * a8;
        a6 = -siniq * a9 + cosiq * a10;
        x1 = a1 * cosomo + a2 * sinomo;
        x2 = a3 * cosomo + a4 * sinomo;
        x3 = -a1 * sinomo + a2 * cosomo;
        x4 = -a3 * sinomo + a4 * cosomo;
        x5 = a5 * sinomo;
        x6 = a6 * sinomo;
        x7 = a5 * cosomo;
        x8 = a6 * cosomo;
        z31 = 12.0 * x1 * x1 - 3. * x3 * x3;
        z32 = 24.0 * x1 * x2 - 6. * x3 * x4;
        z33 = 12.0 * x2 * x2 - 3. * x4 * x4;
        z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * eqsq;
        z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * eqsq;
        z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * eqsq;
        z11 = -6.0 * a1 * a5 + eqsq * (-24. * x1 * x7 - 6. * x3 * x5);
        z12 = -6.0 * (a1 * a6 + a3 * a5) + eqsq * (-24. * (x2 * x7 + x1 * x8) - 6. * (x3 * x6 + x4 * x5));
        z13 = -6.0 * a3 * a6 + eqsq * (-24. * x2 * x8 - 6. * x4 * x6);
        z21 = 6.0 * a2 * a5 + eqsq * (24. * x1 * x5 - 6. * x3 * x7);
        z22 = 6.0 * (a4 * a5 + a2 * a6) + eqsq * (24. * (x2 * x5 + x1 * x6) - 6. * (x4 * x7 + x3 * x8));
        z23 = 6.0 * a4 * a6 + eqsq * (24. * x2 * x6 - 6. * x4 * x8);
        z1 = z1 + z1 + bsq * z31;
        z2 = z2 + z2 + bsq * z32;
        z3 = z3 + z3 + bsq * z33;
        s3 = cc * xnoi;
        s2 = -0.5 * s3 / rteqsq;
        s4 = s3 * rteqsq;
        s1 = -15.0 * Eccentricity() * s4;
        s5 = x1 * x3 + x2 * x4;
        s6 = x2 * x3 + x1 * x4;
        s7 = x2 * x4 - x1 * x3;
        se = s1 * zn * s5;
        si = s2 * zn * (z11 + z13);
        sl = -zn * s3 * (z1 + z3 - 14.0 - 6.0 * eqsq);
        sgh = s4 * zn * (z31 + z33 - 6.0);

        /*
         * replaced
         * sh = -zn * s2 * (z21 + z23
         * with
         * shdq = (-zn * s2 * (z21 + z23)) / siniq
         */
        if (Inclination() < 5.2359877e-2 || Inclination() > Globals::PI() - 5.2359877e-2) {
            shdq = 0.0;
        } else {
            shdq = (-zn * s2 * (z21 + z23)) / siniq;
        }

        d_ee2_ = 2.0 * s1 * s6;
        d_e3_ = 2.0 * s1 * s7;
        d_xi2_ = 2.0 * s2 * z12;
        d_xi3_ = 2.0 * s2 * (z13 - z11);
        d_xl2_ = -2.0 * s3 * z2;
        d_xl3_ = -2.0 * s3 * (z3 - z1);
        d_xl4_ = -2.0 * s3 * (-21.0 - 9.0 * eqsq) * ze;
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
        d_ssg_ = sgh - cosiq * d_ssh_;
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
        zmo = zmol;

    }
    d_sse_ += se;
    d_ssi_ += si;
    d_ssl_ += sl;
    d_ssg_ += sgh - cosiq * shdq;
    d_ssh_ += shdq;

    /*
     * geopotential resonance initialization for 12 hour orbits
     */
    bool resonance_flag = false;
    bool synchronous_flag = false;
    bool initialize_integrator = true;

    if (RecoveredMeanMotion() < 0.0052359877 && RecoveredMeanMotion() > 0.0034906585) {

        if (RecoveredMeanMotion() < 8.26e-3 || RecoveredMeanMotion() > 9.24e-3 || Eccentricity() < 0.5) {
            initialize_integrator = false;
        } else {
            /*
             * geopotential resonance initialization for 12 hour orbits
             */
            resonance_flag = true;

            double eoc = Eccentricity() * eqsq;

            double g211;
            double g310;
            double g322;
            double g410;
            double g422;
            double g520;

            double g201 = -0.306 - (Eccentricity() - 0.64) * 0.440;

            if (Eccentricity() <= 0.65) {
                g211 = 3.616 - 13.247 * Eccentricity() + 16.290 * eqsq;
                g310 = -19.302 + 117.390 * Eccentricity() - 228.419 * eqsq + 156.591 * eoc;
                g322 = -18.9068 + 109.7927 * Eccentricity() - 214.6334 * eqsq + 146.5816 * eoc;
                g410 = -41.122 + 242.694 * Eccentricity() - 471.094 * eqsq + 313.953 * eoc;
                g422 = -146.407 + 841.880 * Eccentricity() - 1629.014 * eqsq + 1083.435 * eoc;
                g520 = -532.114 + 3017.977 * Eccentricity() - 5740 * eqsq + 3708.276 * eoc;
            } else {
                g211 = -72.099 + 331.819 * Eccentricity() - 508.738 * eqsq + 266.724 * eoc;
                g310 = -346.844 + 1582.851 * Eccentricity() - 2415.925 * eqsq + 1246.113 * eoc;
                g322 = -342.585 + 1554.908 * Eccentricity() - 2366.899 * eqsq + 1215.972 * eoc;
                g410 = -1052.797 + 4758.686 * Eccentricity() - 7193.992 * eqsq + 3651.957 * eoc;
                g422 = -3581.69 + 16178.11 * Eccentricity() - 24462.77 * eqsq + 12422.52 * eoc;

                if (Eccentricity() <= 0.715) {
                    g520 = 1464.74 - 4664.75 * Eccentricity() + 3763.64 * eqsq;
                } else {
                    g520 = -5149.66 + 29936.92 * Eccentricity() - 54087.36 * eqsq + 31324.56 * eoc;
                }
            }

            double g533;
            double g521;
            double g532;

            if (Eccentricity() < 0.7) {
                g533 = -919.2277 + 4988.61 * Eccentricity() - 9064.77 * eqsq + 5542.21 * eoc;
                g521 = -822.71072 + 4568.6173 * Eccentricity() - 8491.4146 * eqsq + 5337.524 * eoc;
                g532 = -853.666 + 4690.25 * Eccentricity() - 8624.77 * eqsq + 5341.4 * eoc;
            } else {
                g533 = -37995.78 + 161616.52 * Eccentricity() - 229838.2 * eqsq + 109377.94 * eoc;
                g521 = -51752.104 + 218913.95 * Eccentricity() - 309468.16 * eqsq + 146349.42 * eoc;
                g532 = -40023.88 + 170470.89 * Eccentricity() - 242699.48 * eqsq + 115605.82 * eoc;
            }

            double sini2 = siniq * siniq;
            double f220 = 0.75 * (1.0 + 2.0 * cosiq + cosq2);
            double f221 = 1.5 * sini2;
            double f321 = 1.875 * siniq * (1.0 - 2.0 * cosiq - 3.0 * cosq2);
            double f322 = -1.875 * siniq * (1.0 + 2.0 * cosiq - 3.0 * cosq2);
            double f441 = 35.0 * sini2 * f220;
            double f442 = 39.3750 * sini2 * sini2;
            double f522 = 9.84375 * siniq * (sini2 * (1.0 - 2.0 * cosiq - 5.0 * cosq2)
                    + 0.33333333 * (-2.0 + 4.0 * cosiq + 6.0 * cosq2));
            double f523 = siniq * (4.92187512 * sini2 * (-2.0 - 4.0 * cosiq + 10.0 * cosq2)
                    + 6.56250012 * (1.0 + 2.0 * cosiq - 3.0 * cosq2));
            double f542 = 29.53125 * siniq * (2.0 - 8.0 * cosiq + cosq2 *
                    (-12.0 + 8.0 * cosiq + 10.0 * cosq2));
            double f543 = 29.53125 * siniq * (-2.0 - 8.0 * cosiq + cosq2 *
                    (12.0 + 8.0 * cosiq - 10.0 * cosq2));

            double xno2 = RecoveredMeanMotion() * RecoveredMeanMotion();
            double ainv2 = aqnv * aqnv;

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

            d_xlamo_ = MeanAnomoly() + AscendingNode() + AscendingNode() - gsto_ - gsto_;
            bfact = xlldot + xnodot + xnodot - THDT - THDT;
            bfact = bfact + d_ssl_ + d_ssh_ + d_ssh_;
        }
    } else {
        /*
         * 24h synchronous resonance terms initialization
         */
        resonance_flag = true;
        synchronous_flag = true;

        double g200 = 1.0 + eqsq * (-2.5 + 0.8125 * eqsq);
        double g310 = 1.0 + 2.0 * eqsq;
        double g300 = 1.0 + eqsq * (-6.0 + 6.60937 * eqsq);
        double f220 = 0.75 * (1.0 + cosiq) * (1.0 + cosiq);
        double f311 = 0.9375 * siniq * siniq * (1.0 + 3.0 * cosiq) - 0.75 * (1.0 + cosiq);
        double f330 = 1.0 + cosiq;
        f330 = 1.875 * f330 * f330 * f330;
        d_del1_ = 3.0 * RecoveredMeanMotion() * RecoveredMeanMotion() * aqnv * aqnv;
        d_del2_ = 2.0 * d_del1_ * f220 * g200 * Q22;
        d_del3_ = 3.0 * d_del1_ * f330 * g300 * Q33 * aqnv;
        d_del1_ = d_del1_ * f311 * g310 * Q31 * aqnv;
        d_fasx2_ = 0.13130908;
        d_fasx4_ = 2.8843198;
        d_fasx6_ = 0.37448087;

        d_xlamo_ = MeanAnomoly() + AscendingNode() + ArgumentPerigee() - gsto_;
        bfact = xlldot + xpidot - THDT;
        bfact = bfact + d_ssl_ + d_ssg_ + d_ssh_;
    }

    if (initialize_integrator) {
        d_xfact_ = bfact - RecoveredMeanMotion();
        /*
         * initialize integrator
         */
        d_xli_ = d_xlamo_;
        d_xni_ = RecoveredMeanMotion();
        d_atime_ = 0.0;
        d_stepp_ = 720.0;
        d_stepn_ = -720.0;
        d_step2_ = 259200.0;
    }
}

#if 0
SUBROUTINE DEEP
COMMON / E1 / XMO, XNODEO, OMEGAO, EO, XINCL, XNO, XNDT2O,
        1 XNDD6O, BSTAR, X, Y, Z, XDOT, YDOT, ZDOT, EPOCH, DS50
        COMMON / C1 / CK2, CK4, E6A, QOMS2T, S, TOTHRD,
        1 XJ3, XKE, XKMPER, XMNPDA, AE
        COMMON / C2 / DE2RA, PI, PIO2, TWOPI, X3PIO2
        DOUBLE PRECISION EPOCH, DS50
        DOUBLE PRECISION
        * DAY, PREEP, XNODCE, ATIME, DELT, SAVTSN, STEP2, STEPN, STEPP

        double ZNS = 1.19459E-5;
double C1SS = 2.9864797E-6;
double ZES = 0.01675;
double ZNL = 1.5835218E-4;
double C1L = 4.7968065E-7;
double ZEL = 0.05490;
double ZCOSIS = 0.91744867;
double ZSINI = 0.39785416;
double ZSINGS = -0.98088458;
double ZCOSGS = 0.1945905;
double ZCOSHS = 1.0;
double ZSINHS = 0.0;
double Q22 = 1.7891679E-6;
double Q31 = 2.1460748E-6;
double Q33 = 2.2123015E-7;
double G22 = 5.7686396;
double G32 = 0.95240898;
double G44 = 1.8014998;
double G52 = 1.0508330;
double G54 = 4.4108898;
double ROOT22 = 1.7891679E-6;
double ROOT32 = 3.7393792E-7;
double ROOT44 = 7.3636953E-9;
double ROOT52 = 1.1428639E-7;
double ROOT54 = 2.1765803E-9;
double THDT = 4.3752691E-3;





/*
 * ENTRANCE FOR DEEP SPACE SECULAR EFFECTS
 */
ENTRY DPSEC(XLL, OMGASM, XNODES, EM, XINC, XN, T)
XLL = XLL + SSL*T
OMGASM = OMGASM + SSG*T
XNODES = XNODES + SSH*T
EM = EO + SSE*T
XINC = XINCL + SSI * T
IF(XINC .GE. 0.) GO TO 90
XINC = -XINC
XNODES = XNODES + PI
OMGASM = OMGASM - PI
90 IF(IRESFL .EQ. 0) RETURN
100 IF(ATIME.EQ.0.D0) GO TO 170
IF(T.GE.(0.D0).AND.ATIME.LT.(0.D0)) GO TO 170
IF(T.LT.(0.D0).AND.ATIME.GE.(0.D0)) GO TO 170
105 IF(DABS(T).GE.DABS(ATIME)) GO TO 120
DELT = STEPP
IF(T.GE.0.D0) DELT = STEPN
110 ASSIGN 100 TO IRET
GO TO 160
120 DELT = STEPN
IF(T.GT.0.D0) DELT = STEPP
125 IF(DABS(T - ATIME).LT.STEPP) GO TO 130
ASSIGN 125 TO IRET
GO TO 160
130 FT = T - ATIME
ASSIGN 140 TO IRETN
GO TO 150
140 XN = XNI + XNDOT * FT + XNDDT * FT * FT * 0.5
XL = XLI + XLDOT * FT + XNDOT * FT * FT * 0.5
TEMP = -XNODES + THGR + T*THDT
XLL = XL - OMGASM + TEMP
IF(ISYNFL.EQ.0) XLL = XL + TEMP + TEMP
RETURN
C
C DOT TERMS CALCULATED
C
150 IF(ISYNFL.EQ.0) GO TO 152
XNDOT = DEL1 * SIN(XLI - FASX2) + DEL2 * SIN(2. * (XLI - FASX4))
1 + DEL3 * SIN(3. * (XLI - FASX6))
XNDDT = DEL1 * COS(XLI - FASX2)
* +2. * DEL2 * COS(2. * (XLI - FASX4))
* +3. * DEL3 * COS(3. * (XLI - FASX6))
GO TO 154
152 XOMI = OMEGAQ + OMGDT * ATIME
65
X2OMI = XOMI + XOMI
X2LI = XLI + XLI
XNDOT = D2201 * SIN(X2OMI + XLI - G22)
* +D2211 * SIN(XLI - G22)
* +D3210 * SIN(XOMI + XLI - G32)
* +D3222 * SIN(-XOMI + XLI - G32)
* +D4410 * SIN(X2OMI + X2LI - G44)
* +D4422 * SIN(X2LI - G44)
* +D5220 * SIN(XOMI + XLI - G52)
* +D5232 * SIN(-XOMI + XLI - G52)
* +D5421 * SIN(XOMI + X2LI - G54)
* +D5433 * SIN(-XOMI + X2LI - G54)
XNDDT = D2201 * COS(X2OMI + XLI - G22)
* +D2211 * COS(XLI - G22)
* +D3210 * COS(XOMI + XLI - G32)
* +D3222 * COS(-XOMI + XLI - G32)
* +D5220 * COS(XOMI + XLI - G52)
* +D5232 * COS(-XOMI + XLI - G52)
* +2. * (D4410 * COS(X2OMI + X2LI - G44)
        * +D4422 * COS(X2LI - G44)
        * +D5421 * COS(XOMI + X2LI - G54)
        * +D5433 * COS(-XOMI + X2LI - G54))
154 XLDOT = XNI + XFACT
XNDDT = XNDDT * XLDOT
GO TO IRETN
C
C INTEGRATOR
C
160 ASSIGN 165 TO IRETN
GO TO 150
165 XLI = XLI + XLDOT * DELT + XNDOT*STEP2
XNI = XNI + XNDOT * DELT + XNDDT*STEP2
ATIME = ATIME + DELT
GO TO IRET
C
C EPOCH RESTART
C
170 IF(T.GE.0.D0) GO TO 175
DELT = STEPN
GO TO 180
175 DELT = STEPP
180 ATIME = 0.D0
XNI = XNQ
XLI = XLAMO
GO TO 125




/*
 * ENTRANCES FOR LUNAR-SOLAR PERIODICS
 */
ENTRY DPPER(EM, XINC, OMGASM, XNODES, XLL)
SINIS = SIN(XINC)
COSIS = COS(XINC)
IF(DABS(SAVTSN - T).LT.(30.D0)) GO TO 210
SAVTSN = T
ZM = ZMOS + ZNS * T
205 ZF = ZM + 2. * ZES * SIN(ZM)
SINZF = SIN(ZF)
F2 = .5 * SINZF * SINZF - .25
F3 = -.5 * SINZF * COS(ZF)
SES = SE2 * F2 + SE3*F3
SIS = SI2 * F2 + SI3*F3
SLS = SL2 * F2 + SL3 * F3 + SL4*SINZF
SGHS = SGH2 * F2 + SGH3 * F3 + SGH4*SINZF
SHS = SH2 * F2 + SH3*F3
ZM = ZMOL + ZNL*T
ZF = ZM + 2. * ZEL * SIN(ZM)
SINZF = SIN(ZF)
F2 = .5 * SINZF * SINZF - .25
F3 = -.5 * SINZF * COS(ZF)
SEL = EE2 * F2 + E3*F3
SIL = XI2 * F2 + XI3*F3
SLL = XL2 * F2 + XL3 * F3 + XL4*SINZF
SGHL = XGH2 * F2 + XGH3 * F3 + XGH4*SINZF
SHL = XH2 * F2 + XH3*F3
PE = SES + SEL
PINC = SIS + SIL
PL = SLS + SLL
210 PGH = SGHS + SGHL
PH = SHS + SHL
XINC = XINC + PINC
EM = EM + PE
IF(XQNCL.LT.(.2)) GO TO 220
GO TO 218
C
C APPLY PERIODICS DIRECTLY
C
218 PH = PH / SINIQ
PGH = PGH - COSIQ*PH
OMGASM = OMGASM + PGH
XNODES = XNODES + PH
XLL = XLL + PL
GO TO 230
C
C APPLY PERIODICS WITH LYDDANE MODIFICATION
C
220 SINOK = SIN(XNODES)
67
COSOK = COS(XNODES)
ALFDP = SINIS*SINOK
BETDP = SINIS*COSOK
DALF = PH * COSOK + PINC * COSIS*SINOK
DBET = -PH * SINOK + PINC * COSIS*COSOK
ALFDP = ALFDP + DALF
BETDP = BETDP + DBET
XLS = XLL + OMGASM + COSIS*XNODES
DLS = PL + PGH - PINC * XNODES*SINIS
XLS = XLS + DLS
XNODES = ACTAN(ALFDP, BETDP)
XLL = XLL + PL
OMGASM = XLS - XLL - COS(XINC) * XNODES
230 CONTINUE
RETURN
END
#endif