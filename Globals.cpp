#include "Globals.h"

Globals::Globals(void) {
}

Globals::~Globals(void) {
}

double Globals::Fmod2p(const double& arg) {
    double modu = fmod(arg, TWOPI());
    if (modu < 0.0)
        modu += TWOPI();

    return modu;
}

double Globals::Deg2Rad(const double& deg) {
    return deg * PI() / 180.0;
}

double Globals::Rad2Deg(const double& rad) {
    return rad * 180.0 / PI();
}
