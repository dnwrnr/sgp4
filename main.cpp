#include "Julian.h"
#include "Tle.h"
#include "SGDP4.h"
#include "Globals.h"
#include "Observer.h"
#include "Coord.h"

#include <list>
#include <string>
#include <iomanip>

void RunTest();

void GeneratePassList(const CoordGeodetic& geo, const SGDP4& model, const Julian& date) {
    Observer obs(geo);
    Eci eci;
    model.FindPosition(eci, date);

    CoordTopographic topo = obs.GetLookAngle(eci);

    Julian date_start;
    Julian date_end;

    date_end.AddDay(10.0);

    for (Julian jd = date_start; jd <= date_end; jd.AddMin(1.0)) {

    }
}

int main() {

    Tle tle = Tle("UK-DMC 2                ",
            "1 35683U 09041C   11089.11558659  .00000272  00000-0  54146-4 0  8712",
            "2 35683  98.0762 348.1067 0001434  99.8921 260.2456 14.69414094 89293");

    CoordGeodetic geo(Globals::Deg2Rad(51.360242), Globals::Deg2Rad(0.101473), 0.07);

    SGDP4 sgp4_model;
    sgp4_model.SetTle(tle);
    Julian date;

    GeneratePassList(geo, sgp4_model, date);


    return 0;
}

