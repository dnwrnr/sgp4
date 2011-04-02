#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <cmath>

class Globals {
public:
    Globals(void);
    ~Globals(void);

    static const double F() {
        /*
         * earth flattening
         */
        return 1.0 / 298.26;
    }

    static const double OMEGA_E() {
        /*
         * earth rotation per sideral day
         */
        return 1.00273790934;
    }

    static const double XKMPER() {
        return 6378.135;
    }

    static const double PI() {
        return 3.14159265358979323846264338327950288419716939937510582;
    }

    static const double TWOPI() {
        return 2.0 * 3.14159265358979323846264338327950288419716939937510582;
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

