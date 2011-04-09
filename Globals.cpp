#include "Globals.h"

Globals::Globals(void) {
}

Globals::~Globals(void) {
}

double Globals::Fmod2p(const double arg) {

    double modu = fmod(arg, TWOPI());
    if (modu < 0.0)
        modu += TWOPI();

    return modu;
}

double Globals::Deg2Rad(const double deg) {

    return deg * PI() / 180.0;
}

double Globals::Rad2Deg(const double rad) {

    return rad * 180.0 / PI();
}

double Globals::AcTan(const double sinx, const double cosx) {

    if (cosx == 0.0) {
        if (sinx > 0.0)
            return PI() / 2.0;
        else
            return 3.0 * PI() / 2.0;
    } else {
        if (cosx > 0.0)
            return atan(sinx / cosx);
        else
            return PI() + atan(sinx / cosx);
    }
}