#include "Observer.h"
#include "SGP4.h"
#include "Util.h"

#include <cmath>
#include <iostream>

const int kPASS_TIME_STEP = 180; // time step in seconds, when searching for aos/los

double FindMaxElevation(const CoordGeodetic& user_geo, SGP4& sgp4, const Julian& aos, const Julian& los) {

    Observer obs(user_geo);

    double max_elevation = -99.9;
    /*
     * time step in seconds
     */
    double time_step = 180.0;
    /*
     * still searching for max elevation
     */
    bool searching = true;
    /*
     * fine tune the max elevation value
     */
    bool fine_tune = false;

    Julian current_time = aos;
    while (current_time < los && searching) {

        /*
         * find position
         */
        Eci eci = sgp4.FindPosition(current_time);
        CoordTopographic topo = obs.GetLookAngle(eci);

        /*
         * keep updating max elevation
         */
        if (topo.elevation > max_elevation) {

            max_elevation = topo.elevation;
        } else if (!fine_tune) {
            /*
             * passed max elevation
             * max elevation happened in the last 6 minutes
             * go back and fine tune max elevation value
             */
            current_time.AddSec(-2.0 * time_step);
            /*
             * dont go back before aos
             */
            if (current_time < aos)
                current_time = aos;

            /*
             * 1 second increment
             */
            time_step = 1.0;
            fine_tune = true;

            /*
             * reset elevation
             */
            max_elevation = -99.9;

        } else {
            /*
             * found max elevation
             */

            searching = false;
        }

        if (searching) {
            current_time.AddSec(time_step);
            if (current_time > los)
                current_time = los;
        }
    };

    return max_elevation;
}

Julian FindCrossingPoint(const CoordGeodetic& user_geo, SGP4& sgp4, const Julian& initial_time1, const Julian& initial_time2, bool finding_aos) {

    Observer obs(user_geo);

    bool searching = true;
    unsigned int loop_count = 0;

    Julian time1 = initial_time1;
    Julian time2 = initial_time2;

    Julian middle_time = time1 + (time2 - time1) / 2.0;

    while (searching && loop_count < 25) {

        /*
         * find position
         */
        Eci eci = sgp4.FindPosition(middle_time);
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
        if (((time2.GetDate() - time1.GetDate()) * kSECONDS_PER_DAY) < 0.5) {
            searching = false;
        }

        middle_time = time1 + (time2 - time1) / 2.0;
        loop_count++;
    };

    return middle_time;
}

void AOSLOS(const CoordGeodetic& user_geo, SGP4& sgp4, const Julian& start_time, const Julian& end_time) {

    Observer obs(user_geo);

    Timespan time_step(0, 0, 0, kPASS_TIME_STEP);

    bool first_run = true;
    bool end_of_pass = false;
    Julian previous_time = start_time;

    Julian aos_time;
    Julian los_time;
    bool found_aos = false;
    bool found_los = false;

    Julian current_time = start_time;
    while (current_time < end_time) {

        // find position
        Eci eci = sgp4.FindPosition(current_time);
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
                end_of_pass = true;
                double max_elevation = FindMaxElevation(user_geo, sgp4, aos_time, los_time);
                std::cout << "AOS: " << aos_time << ", LOS: " << los_time << ", Max El: " << Util::RadiansToDegrees(max_elevation) << std::endl;

            }
        }

        first_run = false;
        previous_time = current_time;

        // at end of pass, move time along by 20 minutes
        if (end_of_pass)
            current_time += Timespan(0, 0, 20, 0);
        else
            current_time += time_step;

        // check we dont go past end time
        if (current_time > end_time)
            current_time = end_time;

        end_of_pass = false;
    };

    /*
     * is satellite still above horizon at end of search period
     */
    if (found_aos && !found_los) {

        los_time = end_time;
        double max_elevation = FindMaxElevation(user_geo, sgp4, aos_time, los_time);
        std::cout << "AOS: " << aos_time << ", LOS: " << los_time << ", Max El: " << Util::RadiansToDegrees(max_elevation) << std::endl;
    }
}

int main() {

    CoordGeodetic geo(51.507406923983446, -0.12773752212524414, 0.05);
    Tle tle("GIOVE-B                 ",
            "1 32781U 08020A   11158.03814084  .00000088  00000-0  10000-3 0  4620",
            "2 32781  55.9142 172.9458 0022365 228.3743 131.4697  1.70953903 19437");
    SGP4 sgp4(tle);

    Julian start_date;
    Julian end_date(start_date);
    end_date.AddDay(10.0);

    /*
     * generate 10 day schedule
     */
    AOSLOS(geo, sgp4, start_date, end_date);

    return 0;
}
