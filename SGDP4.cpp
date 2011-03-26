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
    int year = static_cast<int> (tle.GetField(Tle::FLD_EPOCHYEAR));
    if (year < 57)
        year += 2000;
    else
        year += 1900;
    double day = tle.GetField(Tle::FLD_EPOCHDAY);
    Julian jul(year, day);
    epoch_ = jul;


    recovered_semi_major_axis = 0.0;
    recovered_mean_motion_ = 0.0;

    Initialize();
}

void SGDP4::Initialize() {

    cosio_ = 0.0;
    sinio_ = 0.0;
    aodp_ = 0.0;

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

    use_simple_model_ = false;
    use_deep_space_ = false;



    /*
     * recover original mean motion and semimajor axis (aodp)
     * from input elements
     */
    cosio_ = cos(tle_data_0_.xincl);
    sinio_ = sin(tle_data_0_.xincl);
    double theta2 = cosio_ * cosio_;
    x3thm1_ = 3.0 * theta2 - 1.0;
    double eosq = tle_data_0_.eo * tle_data_0_.eo;
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);
    {
        double a1 = pow(Globals::XKE() / tle_data_0_.xno, Globals::TOTHRD());
        double del1 = 1.5 * Globals::CK2() * x3thm1_ / (a1 * a1 * betao * betao2);
        double ao = a1 * (1.0 - del1 * (0.5 * Globals::TOTHRD() + del1 * (1.0 + 134.0 / 81.0 * del1)));
        double delo = 1.5 * Globals::CK2() * x3thm1_ / (ao * ao * betao * betao2);
        /*
         * recovered mean motion
         */
        tle_data_0_.xno = tle_data_0_.xno / (1.0 + delo);
        /*
         * recovered semimajor axis
         */
        aodp_ = ao / (1.0 - delo);
    }

    double rp = aodp_ * (1.0 - tle_data_0_.eo);
    double perigee = (rp - Globals::AE()) * Globals::XKMPER();
    double period = Globals::TWOPI() / tle_data_0_.xno;

    /*
     * for perigee less than 220 kilometers, the simple_model flag is set and
     * the equations are truncated to linear variation in sqrt a and
     * quadratic variation in mean anomly. also, the c3 term, the
     * delta omega term and the delta m term are dropped
     */
    use_deep_space_ = false;
    use_simple_model_ = false;
    if (period >= 225.0) {
        use_deep_space_ = true;
    } else {
        if (rp < (220.0 / Globals::XKMPER() + Globals::AE()))
            use_simple_model_ = true;
    }

    double s4_ = Globals::S();
    double qoms24_ = Globals::QOMS2T();
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

    double pinvsq = 1.0 / (aodp_ * aodp_ * betao2 * betao2);

    double tsi = 1.0 / (aodp_ - s4_);
    eta_ = aodp_ * tle_data_0_.eo * tsi;
    double etasq = eta_ * eta_;
    double eeta = tle_data_0_.eo * eta_;
    double psisq = fabs(1.0 - etasq);
    double coef = qoms24_ * pow(tsi, 4.0);
    coef1_ = coef / pow(psisq, 3.5);
    double c2 = coef1_ * tle_data_0_.xno * (aodp_ * (1.0 + 1.5 * etasq + eeta *
            (4.0 + etasq)) + 0.75 * Globals::CK2() * tsi / psisq *
            x3thm1_ * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    c1_ = tle_data_0_.bstar * c2;
    a3ovk2_ = -Globals::XJ3() / Globals::CK2() * pow(Globals::AE(), 3.0);

    x1mth2_ = 1.0 - theta2;
    c4_ = 2.0 * tle_data_0_.xno * coef1_ * aodp_ * betao2 *
            (eta_ * (2.0 + 0.5 * etasq) + tle_data_0_.eo * (0.5 + 2.0 * etasq) -
            2.0 * Globals::CK2() * tsi / (aodp_ * psisq) *
            (-3.0 * x3thm1_ * (1.0 - 2.0 * eeta + etasq *
            (1.5 - 0.5 * eeta)) + 0.75 * x1mth2_ * (2.0 * etasq - eeta *
            (1.0 + etasq)) * cos(2.0 * tle_data_0_.omega)));
    double theta4 = theta2 * theta2;
    double temp1 = 3.0 * Globals::CK2() * pinvsq * tle_data_0_.xno;
    double temp2 = temp1 * Globals::CK2() * pinvsq;
    double temp3 = 1.25 * Globals::CK4() * pinvsq * pinvsq * tle_data_0_.xno;
    xmdot_ = tle_data_0_.xno + 0.5 * temp1 * betao * x3thm1_ + 0.0625 * temp2 * betao *
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
        if (tle_data_0_.eo > 1.0e-4) {
            c3 = coef * tsi * a3ovk2_ * tle_data_0_.xno * Globals::AE() *
                    sinio_ / tle_data_0_.eo;
        }

        c5_ = 2.0 * coef1_ * aodp_ * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
        omgcof_ = tle_data_0_.bstar * c3 * cos(tle_data_0_.omega);

        xmcof_ = 0.0;
        if (tle_data_0_.eo > 1.0e-4)
            xmcof_ = -Globals::TOTHRD() * coef * tle_data_0_.bstar * Globals::AE() / eeta;

        delmo_ = pow(1.0 + eta_ * (cos(tle_data_0_.xmo)), 3.0);
        sinmo_ = sin(tle_data_0_.xmo);
    }

    if (!use_simple_model_) {
        double c1sq = c1_ * c1_;
        d2_ = 4.0 * aodp_ * tsi * c1sq;
        double temp = d2_ * tsi * c1_ / 3.0;
        d3_ = (17.0 * aodp_ + s4_) * temp;
        d4_ = 0.5 * temp * aodp_ * tsi * (221.0 * aodp_ + 31.0 * s4_) * c1_;
        t3cof_ = d2_ + 2.0 * c1sq;
        t4cof_ = 0.25 * (3.0 * d3_ + c1_ * (12.0 * d2_ + 10.0 * c1sq));
        t5cof_ = 0.2 * (3.0 * d4_ + 12.0 * c1_ * d3_ + 6.0 * d2_ * d2_ + 15.0 *
                c1sq * (2.0 * d2_ + c1sq));
    } else if (use_deep_space_) {
        gsto_ = tle_data_0_.epoch.ToGMST();
        //CALL DPINIT(EOSQ,SINIO,COSIO,BETAO,AODP,THETA2,
        // SING,COSG,BETAO2,XMDOT,OMGDOT,XNODOT,XNODP)
    }

    first_run_ = false;
}

void SGDP4::FindPosition(double tsince) {

    struct TleData tle_data_tsince_;
    memset(&tle_data_tsince_, 0, sizeof (tle_data_tsince_));

    tle_data_tsince_.bstar = tle_data_0_.bstar;
    tle_data_tsince_.eo = tle_data_0_.eo;
    tle_data_tsince_.omega = tle_data_0_.omega;
    tle_data_tsince_.xincl = tle_data_0_.xincl;
    tle_data_tsince_.xmo = tle_data_0_.xmo;
    tle_data_tsince_.xno = tle_data_0_.xno;
    tle_data_tsince_.xnodeo = tle_data_0_.xnodeo;
    tle_data_tsince_.epoch = tle_data_0_.epoch;

    double xl = 0.0;
    double a = 0.0;

    /*
     * update for secular gravity and atmospheric drag
     */
    double xmdf = tle_data_0_.xmo + xmdot_ * tsince;
    double omgadf = tle_data_0_.omega + omgdot_ * tsince;
    double xnoddf = tle_data_tsince_.xnodeo + xnodot_ * tsince;

    double tsq = tsince * tsince;
    double xnode = xnoddf + xnodcf_ * tsq;
    double tempa = 1.0 - c1_ * tsince;
    double tempe = tle_data_tsince_.bstar * c4_ * tsince;
    double templ = t2cof_ * tsq;

    tle_data_tsince_.omega = omgadf;

    if (use_deep_space_) {
        double xn = tle_data_0_.xno;
#if 0
        CALL DPSEC(xmdf, tle_data_tsince_.omega, XNODE, tle_data_tsince_.eo, tle_data_tsince_.xincl, xn, tsince);
#endif
        a = pow(Globals::XKE() / xn, Globals::TOTHRD()) * pow(tempa, 2.0);
        tle_data_tsince_.eo -= tempe;
        double xmam = xmdf + tle_data_0_.xno * templ;
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
        a = aodp_ * pow(tempa, 2.0);
        tle_data_tsince_.eo = tle_data_0_.eo - tempe;
        xl = xmp + tle_data_tsince_.omega + xnode + tle_data_0_.xno * templ;
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
double ZES = .01675;
double ZNL = 1.5835218E-4;
double C1L = 4.7968065E-7;
double ZEL = .05490;
double ZCOSIS = .91744867;
double ZSINI = .39785416;
double ZSINGS = -.98088458;
double ZCOSGS = .1945905;
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
 * ENTRANCE FOR DEEP SPACE INITIALIZATION
 */
ENTRY DPINIT(EQSQ, SINIQ, COSIQ, RTEQSQ, AO, COSQ2, SINOMO, COSOMO,
        1 BSQ, XLLDOT, OMGDT, XNODOT, XNODP)
THGR = THETAG(EPOCH)
EQ = EO
XNQ = XNODP
AQNV = 1. / AO
XQNCL = XINCL
XMAO = XMO
XPIDOT = OMGDT + XNODOT
SINQ = SIN(XNODEO)
COSQ = COS(XNODEO)
OMEGAQ = OMEGAO
* INITIALIZE LUNAR SOLAR TERMS
5 DAY = DS50 + 18261.5D0
IF(DAY.EQ.PREEP) GO TO 10
PREEP = DAY
XNODCE = 4.5236020 - 9.2422029E-4 * DAY
59
STEM = DSIN(XNODCE)
CTEM = DCOS(XNODCE)
ZCOSIL = .91375164 - .03568096 * CTEM
ZSINIL = SQRT(1. - ZCOSIL*ZCOSIL)
ZSINHL = .089683511 * STEM / ZSINIL
ZCOSHL = SQRT(1. - ZSINHL*ZSINHL)
C = 4.7199672 + .22997150 * DAY
GAM = 5.8351514 + .0019443680 * DAY
ZMOL = FMOD2P(C - GAM)
ZX = .39785416 * STEM / ZSINIL
ZY = ZCOSHL*CTEM + 0.91744867 * ZSINHL*STEM
ZX = ACTAN(ZX, ZY)
ZX = GAM + ZX - XNODCE
ZCOSGL = COS(ZX)
ZSINGL = SIN(ZX)
ZMOS = 6.2565837D0 + .017201977D0*DAY
ZMOS = FMOD2P(ZMOS)
* DO SOLAR TERMS
10 LS = 0
SAVTSN = 1.D20
ZCOSG = ZCOSGS
ZSING = ZSINGS
ZCOSI = ZCOSIS
ZSINI = ZSINIS
ZCOSH = COSQ
ZSINH = SINQ
CC = C1SS
ZN = ZNS
ZE = ZES
ZMO = ZMOS
XNOI = 1. / XNQ
ASSIGN 30 TO LS
20 A1 = ZCOSG*ZCOSH + ZSING*ZCOSI*ZSINH
A3 = -ZSING*ZCOSH + ZCOSG*ZCOSI*ZSINH
A7 = -ZCOSG*ZSINH + ZSING*ZCOSI*ZCOSH
A8 = ZSING*ZSINI
A9 = ZSING*ZSINH + ZCOSG*ZCOSI*ZCOSH
A10 = ZCOSG*ZSINI
A2 = COSIQ*A7 + SINIQ*A8
A4 = COSIQ*A9 + SINIQ*A10
A5 = -SINIQ*A7 + COSIQ*A8
A6 = -SINIQ*A9 + COSIQ*A10
C
X1 = A1*COSOMO + A2*SINOMO
X2 = A3*COSOMO + A4*SINOMO
X3 = -A1*SINOMO + A2*COSOMO
60
X4 = -A3*SINOMO + A4*COSOMO
X5 = A5*SINOMO
X6 = A6*SINOMO
X7 = A5*COSOMO
X8 = A6*COSOMO
C
Z31 = 12. * X1*X1 - 3. * X3*X3
Z32 = 24. * X1*X2 - 6. * X3*X4
Z33 = 12. * X2*X2 - 3. * X4*X4
Z1 = 3. * (A1*A1 + A2*A2) + Z31*EQSQ
Z2 = 6. * (A1*A3 + A2*A4) + Z32*EQSQ
Z3 = 3. * (A3*A3 + A4*A4) + Z33*EQSQ
Z11 = -6. * A1*A5 + EQSQ *(-24. * X1*X7 - 6. * X3*X5)
Z12 = -6. * (A1*A6 + A3*A5) + EQSQ *(-24. * (X2*X7 + X1*X8) - 6. * (X3*X6 + X4*X5))
Z13 = -6. * A3*A6 + EQSQ *(-24. * X2*X8 - 6. * X4*X6)
Z21 = 6. * A2*A5 + EQSQ *(24. * X1*X5 - 6. * X3*X7)
Z22 = 6. * (A4*A5 + A2*A6) + EQSQ *(24. * (X2*X5 + X1*X6) - 6. * (X4*X7 + X3*X8))
Z23 = 6. * A4*A6 + EQSQ *(24. * X2*X6 - 6. * X4*X8)
Z1 = Z1 + Z1 + BSQ*Z31
Z2 = Z2 + Z2 + BSQ*Z32
Z3 = Z3 + Z3 + BSQ*Z33
S3 = CC*XNOI
S2 = -.5 * S3 / RTEQSQ
S4 = S3*RTEQSQ
S1 = -15. * EQ*S4
S5 = X1*X3 + X2*X4
S6 = X2*X3 + X1*X4
S7 = X2*X4 - X1*X3
SE = S1*ZN*S5
SI = S2*ZN*(Z11 + Z13)
SL = -ZN*S3*(Z1 + Z3 - 14. - 6. * EQSQ)
SGH = S4*ZN*(Z31 + Z33 - 6.)
SH = -ZN*S2*(Z21 + Z23)
IF(XQNCL.LT.5.2359877E-2) SH = 0.0
EE2 = 2. * S1*S6
E3 = 2. * S1*S7
XI2 = 2. * S2*Z12
XI3 = 2. * S2*(Z13 - Z11)
XL2 = -2. * S3*Z2
XL3 = -2. * S3*(Z3 - Z1)
XL4 = -2. * S3*(-21. - 9. * EQSQ) * ZE
XGH2 = 2. * S4*Z32
XGH3 = 2. * S4*(Z33 - Z31)
XGH4 = -18. * S4*ZE
XH2 = -2. * S2*Z22
XH3 = -2. * S2*(Z23 - Z21)
GO TO LS
61
* DO LUNAR TERMS
30 SSE = SE
SSI = SI
SSL = SL
SSH = SH / SINIQ
SSG = SGH - COSIQ*SSH
SE2 = EE2
SI2 = XI2
SL2 = XL2
SGH2 = XGH2
SH2 = XH2
SE3 = E3
SI3 = XI3
SL3 = XL3
SGH3 = XGH3
SH3 = XH3
SL4 = XL4
SGH4 = XGH4
LS = 1
ZCOSG = ZCOSGL
ZSING = ZSINGL
ZCOSI = ZCOSIL
ZSINI = ZSINIL
ZCOSH = ZCOSHL*COSQ + ZSINHL*SINQ
ZSINH = SINQ*ZCOSHL - COSQ*ZSINHL
ZN = ZNL
CC = C1L
ZE = ZEL
ZMO = ZMOL
ASSIGN 40 TO LS
GO TO 20
40 SSE = SSE + SE
SSI = SSI + SI
SSL = SSL + SL
SSG = SSG + SGH - COSIQ / SINIQ*SH
SSH = SSH + SH / SINIQ
* GEOPOTENTIAL RESONANCE INITIALIZATION FOR 12 HOUR ORBITS
IRESFL = 0
ISYNFL = 0
IF(XNQ.LT.(.0052359877).AND.XNQ.GT.(.0034906585)) GO TO 70
IF(XNQ.LT.(8.26E-3) .OR. XNQ.GT.(9.24E-3)) RETURN
IF(EQ.LT.0.5) RETURN
IRESFL = 1
EOC = EQ*EQSQ
G201 = -.306 - (EQ - .64)*.440
62
IF(EQ.GT.(.65)) GO TO 45
G211 = 3.616 - 13.247 * EQ + 16.290 * EQSQ
G310 = -19.302 + 117.390 * EQ - 228.419 * EQSQ + 156.591 * EOC
G322 = -18.9068 + 109.7927 * EQ - 214.6334 * EQSQ + 146.5816 * EOC
G410 = -41.122 + 242.694 * EQ - 471.094 * EQSQ + 313.953 * EOC
G422 = -146.407 + 841.880 * EQ - 1629.014 * EQSQ + 1083.435 * EOC
G520 = -532.114 + 3017.977 * EQ - 5740 * EQSQ + 3708.276 * EOC
GO TO 55
45 G211 = -72.099 + 331.819 * EQ - 508.738 * EQSQ + 266.724 * EOC
G310 = -346.844 + 1582.851 * EQ - 2415.925 * EQSQ + 1246.113 * EOC
G322 = -342.585 + 1554.908 * EQ - 2366.899 * EQSQ + 1215.972 * EOC
G410 = -1052.797 + 4758.686 * EQ - 7193.992 * EQSQ + 3651.957 * EOC
G422 = -3581.69 + 16178.11 * EQ - 24462.77 * EQSQ + 12422.52 * EOC
IF(EQ.GT.(.715)) GO TO 50
G520 = 1464.74 - 4664.75 * EQ + 3763.64 * EQSQ
GO TO 55
50 G520 = -5149.66 + 29936.92 * EQ - 54087.36 * EQSQ + 31324.56 * EOC
55 IF(EQ.GE.(.7)) GO TO 60
G533 = -919.2277 + 4988.61 * EQ - 9064.77 * EQSQ + 5542.21 * EOC
G521 = -822.71072 + 4568.6173 * EQ - 8491.4146 * EQSQ + 5337.524 * EOC
G532 = -853.666 + 4690.25 * EQ - 8624.77 * EQSQ + 5341.4 * EOC
GO TO 65
60 G533 = -37995.78 + 161616.52 * EQ - 229838.2 * EQSQ + 109377.94 * EOC
G521 = -51752.104 + 218913.95 * EQ - 309468.16 * EQSQ + 146349.42 * EOC
G532 = -40023.88 + 170470.89 * EQ - 242699.48 * EQSQ + 115605.82 * EOC
65 SINI2 = SINIQ*SINIQ
F220 = .75 * (1. + 2. * COSIQ + COSQ2)
F221 = 1.5 * SINI2
F321 = 1.875 * SINIQ*(1. - 2. * COSIQ - 3. * COSQ2)
F322 = -1.875 * SINIQ*(1. + 2. * COSIQ - 3. * COSQ2)
F441 = 35. * SINI2*F220
F442 = 39.3750 * SINI2*SINI2
F522 = 9.84375 * SINIQ*(SINI2*(1. - 2. * COSIQ - 5. * COSQ2)
        1 + .33333333 * (-2. + 4. * COSIQ + 6. * COSQ2))
F523 = SINIQ*(4.92187512 * SINI2*(-2. - 4. * COSIQ + 10. * COSQ2)
        * +6.56250012 * (1. + 2. * COSIQ - 3. * COSQ2))
F542 = 29.53125 * SINIQ*(2. - 8. * COSIQ + COSQ2*(-12. + 8. * COSIQ
        * +10. * COSQ2))
F543 = 29.53125 * SINIQ*(-2. - 8. * COSIQ + COSQ2*(12. + 8. * COSIQ - 10. * COSQ2))
XNO2 = XNQ*XNQ
AINV2 = AQNV*AQNV
TEMP1 = 3. * XNO2*AINV2
TEMP = TEMP1*ROOT22
D2201 = TEMP*F220*G201
D2211 = TEMP*F221*G211
TEMP1 = TEMP1*AQNV
TEMP = TEMP1*ROOT32
D3210 = TEMP*F321*G310
63
D3222 = TEMP*F322*G322
TEMP1 = TEMP1*AQNV
TEMP = 2. * TEMP1*ROOT44
D4410 = TEMP*F441*G410
D4422 = TEMP*F442*G422
TEMP1 = TEMP1*AQNV
TEMP = TEMP1*ROOT52
D5220 = TEMP*F522*G520
D5232 = TEMP*F523*G532
TEMP = 2. * TEMP1*ROOT54
D5421 = TEMP*F542*G521
D5433 = TEMP*F543*G533
XLAMO = XMAO + XNODEO + XNODEO - THGR - THGR
BFACT = XLLDOT + XNODOT + XNODOT - THDT - THDT
BFACT = BFACT + SSL + SSH + SSH
GO TO 80
* SYNCHRONOUS RESONANCE TERMS INITIALIZATION
70 IRESFL = 1
ISYNFL = 1
G200 = 1.0 + EQSQ*(-2.5 + .8125 * EQSQ)
G310 = 1.0 + 2.0 * EQSQ
G300 = 1.0 + EQSQ*(-6.0 + 6.60937 * EQSQ)
F220 = .75 * (1. + COSIQ)*(1. + COSIQ)
F311 = .9375 * SINIQ*SINIQ*(1. + 3. * COSIQ) - .75 * (1. + COSIQ)
F330 = 1. + COSIQ
F330 = 1.875 * F330*F330*F330
DEL1 = 3. * XNQ*XNQ*AQNV*AQNV
DEL2 = 2. * DEL1*F220*G200*Q22
DEL3 = 3. * DEL1*F330*G300*Q33*AQNV
DEL1 = DEL1*F311*G310*Q31*AQNV
FASX2 = .13130908
FASX4 = 2.8843198
FASX6 = .37448087
XLAMO = XMAO + XNODEO + OMEGAO - THGR
BFACT = XLLDOT + XPIDOT - THDT
BFACT = BFACT + SSL + SSG + SSH
80 XFACT = BFACT - XNQ
C
C INITIALIZE INTEGRATOR
C
XLI = XLAMO
XNI = XNQ
ATIME = 0.D0
STEPP = 720.D0
STEPN = -720.D0
STEP2 = 259200.D0
64
RETURN





/*
 * ENTRANCE FOR DEEP SPACE SECULAR EFFECTS
 */
ENTRY DPSEC(XLL, OMGASM, XNODES, EM, XINC, XN, T)
XLL = XLL + SSL*T
OMGASM = OMGASM + SSG*T
XNODES = XNODES + SSH*T
EM = EO + SSE*T
XINC = XINCL + SSI*T
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
140 XN = XNI + XNDOT*FT + XNDDT*FT*FT * 0.5
XL = XLI + XLDOT*FT + XNDOT*FT*FT * 0.5
TEMP = -XNODES + THGR + T*THDT
XLL = XL - OMGASM + TEMP
IF(ISYNFL.EQ.0) XLL = XL + TEMP + TEMP
RETURN
C
C DOT TERMS CALCULATED
C
150 IF(ISYNFL.EQ.0) GO TO 152
XNDOT = DEL1*SIN(XLI - FASX2) + DEL2*SIN(2. * (XLI - FASX4))
1 + DEL3*SIN(3. * (XLI - FASX6))
XNDDT = DEL1*COS(XLI - FASX2)
* +2. * DEL2*COS(2. * (XLI - FASX4))
* +3. * DEL3*COS(3. * (XLI - FASX6))
GO TO 154
152 XOMI = OMEGAQ + OMGDT*ATIME
65
X2OMI = XOMI + XOMI
X2LI = XLI + XLI
XNDOT = D2201*SIN(X2OMI + XLI - G22)
* +D2211*SIN(XLI - G22)
* +D3210*SIN(XOMI + XLI - G32)
* +D3222*SIN(-XOMI + XLI - G32)
* +D4410*SIN(X2OMI + X2LI - G44)
* +D4422*SIN(X2LI - G44)
* +D5220*SIN(XOMI + XLI - G52)
* +D5232*SIN(-XOMI + XLI - G52)
* +D5421*SIN(XOMI + X2LI - G54)
* +D5433*SIN(-XOMI + X2LI - G54)
XNDDT = D2201*COS(X2OMI + XLI - G22)
* +D2211*COS(XLI - G22)
* +D3210*COS(XOMI + XLI - G32)
* +D3222*COS(-XOMI + XLI - G32)
* +D5220*COS(XOMI + XLI - G52)
* +D5232*COS(-XOMI + XLI - G52)
* +2. * (D4410*COS(X2OMI + X2LI - G44)
        * +D4422*COS(X2LI - G44)
        * +D5421*COS(XOMI + X2LI - G54)
        * +D5433*COS(-XOMI + X2LI - G54))
154 XLDOT = XNI + XFACT
XNDDT = XNDDT*XLDOT
GO TO IRETN
C
C INTEGRATOR
C
160 ASSIGN 165 TO IRETN
GO TO 150
165 XLI = XLI + XLDOT*DELT + XNDOT*STEP2
XNI = XNI + XNDOT*DELT + XNDDT*STEP2
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
ZM = ZMOS + ZNS*T
205 ZF = ZM + 2. * ZES*SIN(ZM)
SINZF = SIN(ZF)
F2 = .5 * SINZF*SINZF - .25
F3 = -.5 * SINZF*COS(ZF)
SES = SE2*F2 + SE3*F3
SIS = SI2*F2 + SI3*F3
SLS = SL2*F2 + SL3*F3 + SL4*SINZF
SGHS = SGH2*F2 + SGH3*F3 + SGH4*SINZF
SHS = SH2*F2 + SH3*F3
ZM = ZMOL + ZNL*T
ZF = ZM + 2. * ZEL*SIN(ZM)
SINZF = SIN(ZF)
F2 = .5 * SINZF*SINZF - .25
F3 = -.5 * SINZF*COS(ZF)
SEL = EE2*F2 + E3*F3
SIL = XI2*F2 + XI3*F3
SLL = XL2*F2 + XL3*F3 + XL4*SINZF
SGHL = XGH2*F2 + XGH3*F3 + XGH4*SINZF
SHL = XH2*F2 + XH3*F3
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
DALF = PH*COSOK + PINC*COSIS*SINOK
DBET = -PH*SINOK + PINC*COSIS*COSOK
ALFDP = ALFDP + DALF
BETDP = BETDP + DBET
XLS = XLL + OMGASM + COSIS*XNODES
DLS = PL + PGH - PINC*XNODES*SINIS
XLS = XLS + DLS
XNODES = ACTAN(ALFDP, BETDP)
XLL = XLL + PL
OMGASM = XLS - XLL - COS(XINC) * XNODES
230 CONTINUE
RETURN
END
#endif