#include "Observer.h"
#include "SGP4.h"

#include <cmath>
#include <iostream>

Julian FindCrossingPoint(const CoordGeodetic& user_geo, SGP4& sgp4, const Julian& initial_time1, const Julian& initial_time2, bool finding_aos) {

    Observer obs(user_geo);
    Eci eci;

    bool searching = true;
    unsigned int loop_count = 0;

    Julian time1 = initial_time1;
    Julian time2 = initial_time2;

    Julian middle_time = time1 + (time2.GetDate() - time1.GetDate()) / 2.0;

    while (searching && loop_count < 25) {

        /*
         * find position
         */
        sgp4.FindPosition(&eci, middle_time);
        CoordTopographic topo = obs.GetLookAngle(eci);

        if (topo.elevation > 0.0) {

            if (finding_aos) {
                time2 = middle_time;
            } else {
                time1 = middle_time;
            }
        } else {

            if (finding_aos) {
                time1 = middle_time;
            } else {
                time2 = middle_time;
            }
        }

        /*
         * when two times are within a second, stop
         */
        if (((time2.GetDate() - time1.GetDate()) * kSECONDS_PER_DAY) < 1.0) {
            searching = false;
        }

        middle_time = time1 + (time2.GetDate() - time1.GetDate()) / 2.0;
        loop_count++;
    };

    return middle_time;
}

void AOSLOS(const CoordGeodetic& user_geo, SGP4& sgp4, const Julian& start_time, const Julian& end_time) {

    Observer obs(user_geo);
    Eci eci;

    bool first_run = true;
    Julian previous_time = start_time;

    Julian aos_time;
    Julian los_time;
    bool found_aos = false;
    bool found_los = false;

    for (Julian current_time = start_time; current_time <= end_time; current_time.AddMin(1.0)) {

        sgp4.FindPosition(&eci, current_time);
        CoordTopographic topo = obs.GetLookAngle(eci);

        if (topo.elevation > 0.0) {
            /*
             * satellite is above horizon
             * and its the first run, so just use start_time
             */
            if (first_run) {

                aos_time = start_time;
                found_aos = true;
                found_los = false;

                /*
                 * not the first iteration and satellite has
                 * gone above horizon, so find AOS
                 */
            } else if (!found_aos) {

                aos_time = FindCrossingPoint(user_geo, sgp4, previous_time, current_time, true);
                found_aos = true;
                found_los = false;
            }
        } else {
            /*
             * if satellite now below horizon and have previously
             * found AOS, find LOS
             */
            if (found_aos) {

                los_time = FindCrossingPoint(user_geo, sgp4, previous_time, current_time, false);
                found_aos = false;
                found_los = false;
                std::cout << "AOS: " << aos_time << ", LOS: " << los_time << std::endl;
            }
        }

        first_run = false;
        previous_time = current_time;
    }

    /*
     * is satellite still above horizon at end of search period
     */
    if (found_aos && !found_los) {

        los_time = end_time;
        std::cout << "AOS: " << aos_time << ", LOS: " << los_time << std::endl;
    }
}

int main() {

    CoordGeodetic geo(51.37322, 0.089607, 0.05);
    Tle tle("STS 134                 ",
            "1 37577U 11020A   11146.50425755  .00016213  00000-0  10821-3 0   213",
            "2 37577  51.6484 271.9502 0003400 321.9646 127.1023 15.75523823  1551");
    SGP4 sgp4(tle);

    Julian start_date;
    Julian end_date(start_date);
    end_date.AddDay(10.0);

    /*
     * generate 10 day schedule
     */
    AOSLOS(geo, sgp4, start_date, end_date);

    /*
     * convert time to whole seconds
     */
    //if(finding_aos)
    //    middle_time = floor(middle_time.GetDate() * kSECONDS_PER_DAY) / kSECONDS_PER_DAY;
    //else
    //    middle_time = ceil(middle_time.GetDate() * kSECONDS_PER_DAY) / kSECONDS_PER_DAY;

    return 0;
#if 0
    Eci eci;

    while (1) {
        Julian now;
        sgp4.FindPosition(&eci, now);
        CoordTopographic topo = obs.GetLookAngle(eci);
        CoordGeodetic geo = eci.ToGeodetic();
        std::cout << now << " ";
        std::cout << topo << " ";
        std::cout << geo << std::endl;
    };

    return 0;
#endif
}
