#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <cmath>

class Globals {
public:
    Globals(void);
    ~Globals(void);

    static const double Q0() {
        return 120.0;
    }

    static const double S0() {
        return 78.0;
    }

    static const double CK2() {
        return 0.5 * XJ2() * AE() * AE();
    }

    static const double CK4() {
        return -0.375 * XJ4() * AE() * AE() * AE() * AE();
    }

    static const double TOTHRD() {
        return 2.0 / 3.0;
    }

    static const double XKE() {
        return 60.0 / sqrt(XKMPER() * XKMPER() * XKMPER() / MU());
    }

    static const double S() {
        return AE() * (1.0 + S0() / XKMPER());
    }

    static const double AE() {
        return 1.0;
    }

    static const double XKMPER() {
        return 6378.137;
    }

    static const double QOMS2T() {
        return pow((Q0() - S0()) * AE() / XKMPER(), 4.0);
    }

    static const double XJ2() {
        return 1.08262998905e-3;
    }

    static const double XJ3() {
        return -2.53215306e-6;
    }

    static const double XJ4() {
        return -1.61098761e-6;
    }

    static const double MU() {
        return 398600.8;
    }

    static const double J3OJ2() {
        return XJ3() / XJ2();
    }

    static const double RADS_PER_DEG() {
        return -0.253881E-5;
    }

    static const double PI() {
        return 3.141592653589793238462643383279;
    }

    static const double TWOPI() {
        return 2.0 * PI();
    }

    static const double SEC_PER_DAY() {
        return 86400.0;
    }

    static const double MIN_PER_DAY() {
        return 1440.0;
    }

    static const double HR_PER_DAY() {
        return 24.0;
    }

    static const double OMEGA_E() {
        return 1.00273790934;
    }

    static const double F() {
        return 1.0 / 298.26;
    }

    static const double EPOCH_JAN1_00H_1900() {
        // Jan 1.0 1900 = Jan 1 1900 00h UTC
        return 2415019.5;
    }

    static const double EPOCH_JAN1_12H_1900() {
        // Jan 1.5 1900 = Jan 1 1900 12h UTC
        return 2415020.0;
    }

    static const double EPOCH_JAN1_12H_2000() {
        // Jan 1.5 2000 = Jan 1 2000 12h UTC
        return 2451545.0;
    }

    static double Fmod2p(const double& arg);
    static double Deg2Rad(const double& deg);
    static double Rad2Deg(const double& rad);
};

#endif

