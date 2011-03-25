#include "SGDP4.h"

#include "SatelliteException.h"

#include <math.h>

SGDP4::SGDP4(void) {
    first_run_ = true;
}

SGDP4::~SGDP4(void) {
}

void SGDP4::Initialize(const Tle& tle) {

    cosio_ = 0.0;
    sinio_ = 0.0;
    xnodp_ = 0.0;
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
     * extract and format tle data
     */
    tle_data_0_.bstar = tle.GetField(Tle::FLD_BSTAR);
    tle_data_0_.eo = tle.GetField(Tle::FLD_E);
    tle_data_0_.omega = tle.GetField(Tle::FLD_ARGPER, Tle::U_RAD);
    tle_data_0_.xincl = tle.GetField(Tle::FLD_I, Tle::U_RAD);
    tle_data_0_.xmo = tle.GetField(Tle::FLD_M, Tle::U_RAD);
    tle_data_0_.xno = tle.GetField(Tle::FLD_MMOTION) / (1440.0 / Globals::TWOPI());
    tle_data_0_.xnodeo = tle.GetField(Tle::FLD_RAAN, Tle::U_RAD);

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
    tle_data_0_.epoch = jul;

    /*
     * recover original mean motion (xnodp) and semimajor axis (aodp)
     * from input elements
     */
    double a1 = pow(Globals::XKE() / tle_data_0_.xno, Globals::TOTHRD());
    cosio_ = cos(tle_data_0_.xincl);
    sinio_ = sin(tle_data_0_.xincl);
    double theta2 = cosio_ * cosio_;
    x3thm1_ = 3.0 * theta2 - 1.0;
    double eosq = tle_data_0_.eo * tle_data_0_.eo;
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);
    double del1 = 1.5 * Globals::CK2() * x3thm1_ / (a1 * a1 * betao * betao2);
    double ao = a1 * (1.0 - del1 * (0.5 * Globals::TOTHRD() + del1 * (1.0 + 134.0 / 81.0 * del1)));
    double delo = 1.5 * Globals::CK2() * x3thm1_ / (ao * ao * betao * betao2);
    /*
     * recovered mean motion
     */
    xnodp_ = tle_data_0_.xno / (1.0 + delo);
    /*
     * recovered semimajor axis
     */
    aodp_ = ao / (1.0 - delo);

    gsto_ = tle_data_0_.epoch.ToGMST();

    double rp = aodp_ * (1.0 - tle_data_0_.eo);
    double perigee = (aodp_ * (1.0 - tle_data_0_.eo) - Globals::AE()) * Globals::XKMPER();
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

    double pinvsq = 1.0 / (aodp_ * aodp_ * betao2 * betao2);

    double tsi = 1.0 / (aodp_ - s4_);
    eta_ = aodp_ * tle_data_0_.eo * tsi;
    double etasq = eta_ * eta_;
    double eeta = tle_data_0_.eo * eta_;
    double psisq = fabs(1.0 - etasq);
    double coef = qoms24_ * pow(tsi, 4.0);
    coef1_ = coef / pow(psisq, 3.5);
    double c2 = coef1_ * xnodp_ * (aodp_ * (1.0 + 1.5 * etasq + eeta *
            (4.0 + etasq)) + 0.75 * Globals::CK2() * tsi / psisq *
            x3thm1_ * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    c1_ = tle_data_0_.bstar * c2;
    a3ovk2_ = -Globals::XJ3() / Globals::CK2() * pow(Globals::AE(), 3.0);

    x1mth2_ = 1.0 - theta2;
    c4_ = 2.0 * xnodp_ * coef1_ * aodp_ * betao2 *
            (eta_ * (2.0 + 0.5 * etasq) + tle_data_0_.eo * (0.5 + 2.0 * etasq) -
            2.0 * Globals::CK2() * tsi / (aodp_ * psisq) *
            (-3.0 * x3thm1_ * (1.0 - 2.0 * eeta + etasq *
            (1.5 - 0.5 * eeta)) + 0.75 * x1mth2_ * (2.0 * etasq - eeta *
            (1.0 + etasq)) * cos(2.0 * tle_data_0_.omega)));
    double theta4 = theta2 * theta2;
    double temp1 = 3.0 * Globals::CK2() * pinvsq * xnodp_;
    double temp2 = temp1 * Globals::CK2() * pinvsq;
    double temp3 = 1.25 * Globals::CK4() * pinvsq * pinvsq * xnodp_;
    xmdot_ = xnodp_ + 0.5 * temp1 * betao * x3thm1_ + 0.0625 * temp2 * betao *
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

        // check result (different to telsko)
        double c3 = 0.0;
        if (tle_data_0_.eo > 1.0e-4) {
            c3 = coef * tsi * a3ovk2_ * xnodp_ * Globals::AE() *
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
        double xn = xnodp_;
#if 0
        CALL DPSEC(xmdf, tle_data_tsince_.omega, XNODE, tle_data_tsince_.eo, tle_data_tsince_.xincl, xn, tsince);
#endif
        a = pow(Globals::XKE() / xn, Globals::TOTHRD()) * pow(tempa, 2.0);
        tle_data_tsince_.eo -= tempe;
        double xmam = xmdf + xnodp_ * templ;
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
        if (use_simple_model_) {
            double delomg = omgcof_ * tsince;
            double delm = xmcof_ * (pow(1.0 + eta_ * cos(xmdf), 3.0) - delmo_);
            double temp = delomg + delm;
            xmp = xmdf + temp;
            tle_data_tsince_.omega -= temp;
            double tcube = tsq * tsince;
            double tfour = tsince * tcube;
            tempa -= d2_ * tsq - d3_ * tcube - d4_ * tfour;
            tempe += tle_data_tsince_.bstar * c5_ * (sin(xmp) - sinmo_);
            templ += t3cof_ * tcube + tfour * (t4cof_ + tsince * t5cof_);
        }
        a = aodp_ * pow(tempa, 2.0);
        tle_data_tsince_.eo = tle_data_0_.eo - tempe;
        xl = xmp + tle_data_tsince_.omega + xnode + xnodp_ * templ;
    }

    double beta = sqrt(1.0 - tle_data_tsince_.eo * tle_data_tsince_.eo);
    double xn = Globals::XKE() / pow(a, 1.5);
    /*
     * long period periodics
     */
    double axn = tle_data_tsince_.eo * cos(tle_data_tsince_.omega);
    double temp = 1.0 / (a * beta * beta);
    double xll = temp * xlcof_ * axn;
    double aynl = temp * aycof_;
    double xlt = xl + xll;
    double ayn = tle_data_tsince_.eo * sin(tle_data_tsince_.omega) + aynl;
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
    double elsq = axn * axn + ayn * ayn;
    temp = 1.0 - elsq;
    double pl = a * temp;
    double r = a * (1.0 - ecose);
    double temp1 = 1.0 / r;
    double rdot = Globals::XKE() * sqrt(a) * esine * temp1;
    double rfdot = Globals::XKE() * sqrt(pl) * temp1;
    double temp2 = a * temp1;
    double betal = sqrt(temp);
    double temp3 = 1.0 / (1.0 + betal);
    double cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
    double sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
    double u = atan2(sinu, cosu);
    double sin2u = 2.0 * sinu * cosu;
    double cos2u = 2.0 * cosu * cosu - 1.0;

    temp = 1.0 / pl;
    temp1 = Globals::CK2() * temp;
    temp2 = temp1 * temp;
    /*
     * update for short periodics
     */
    double rk = r * (1.0 - 1.5 * temp2 * betal * x3thm1_) + 0.5 * temp1 * x1mth2_ * cos2u;
    double uk = u - 0.25 * temp2 * x7thm1_ * sin2u;
    double xnodek = xnode + 1.5 * temp2 * cosio_ * sin2u;

    double xinck = tle_data_tsince_.xincl + 1.5 * temp2 * cosio_ * sinio_ * cos2u;
    double rdotk = rdot - xn * temp1 * x1mth2_ * sin2u;
    double rfdotk = rfdot + xn * temp1 * (x1mth2_ * cos2u + 1.5 * x3thm1_);
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
    double x = rk * ux;
    double y = rk * uy;
    double z = rk * uz;
    double xdot = rdotk * ux + rfdotk * vx;
    double ydot = rdotk * uy + rfdotk * vy;
    double zdot = rdotk * uz + rfdotk * vz;
}

#if 0
* ENTRANCE FOR DEEP SPACE SECULAR EFFECTS
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




c
c entrances for lunar - solar periodics
    c
        66
        c
        entry dpper(em, xinc, omgasm, xnodes, xll)
    sinis = sin(xinc)
    cosis = cos(xinc)
    if (dabs(savtsn - t).lt.(30.d0)) go to 210
            savtsn = t
            zm = zmos + zns * t
            205 zf = zm + 2. * zes * sin(zm)
        sinzf = sin(zf)
        f2 = .5 * sinzf * sinzf - .25
            f3 = -.5 * sinzf * cos(zf)
        ses = se2 * f2 + se3 * f3
            sis = si2 * f2 + si3 * f3
            sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf
            sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf
            shs = sh2 * f2 + sh3 * f3
            zm = zmol + znl * t
            zf = zm + 2. * zel * sin(zm)
        sinzf = sin(zf)
        f2 = .5 * sinzf * sinzf - .25
            f3 = -.5 * sinzf * cos(zf)
        sel = ee2 * f2 + e3 * f3
            sil = xi2 * f2 + xi3 * f3
            sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf
            sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf
            shl = xh2 * f2 + xh3 * f3
            pe = ses + sel
            pinc = sis + sil
            pl = sls + sll
            210 pgh = sghs + sghl
            ph = shs + shl
            xinc = xinc + pinc
            em = em + pe
        if (xqncl.lt.(.2)) go to 220
                go to 218
                c
                c apply periodics directly
                c
                218 ph = ph / siniq
                pgh = pgh - cosiq * ph
                omgasm = omgasm + pgh
                xnodes = xnodes + ph
                xll = xll + pl
                go to 230
                c
                c apply periodics with lyddane modification
                c
                220 sinok = sin(xnodes)
            67
            cosok = cos(xnodes)
            alfdp = sinis * sinok
                betdp = sinis * cosok
                dalf = ph * cosok + pinc * cosis * sinok
                dbet = -ph * sinok + pinc * cosis * cosok
                alfdp = alfdp + dalf
                betdp = betdp + dbet
                xls = xll + omgasm + cosis * xnodes
                dls = pl + pgh - pinc * xnodes * sinis
                xls = xls + dls
                xnodes = actan(alfdp, betdp)
            xll = xll + pl
                omgasm = xls - xll - cos(xinc) * xnodes
                230 continue
        }
#endif