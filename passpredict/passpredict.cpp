#include <Observer.h>
#include <SGP4.h>
#include <Util.h>
#include <CoordTopographic.h>
#include <CoordGeodetic.h>

#include <cmath>
#include <iostream>
#include <list>

struct PassDetails
{
    DateTime aos;
    DateTime los;
    double max_elevation;
};

double FindMaxElevation(
        const CoordGeodetic& user_geo,
        SGP4& sgp4,
        const DateTime& aos,
        const DateTime& los)
{
    Observer obs(user_geo);

    double max_elevation = 0.0;

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

    DateTime current_time(aos);
    while (current_time < los && searching)
    {
        /*
         * find position
         */
        Eci eci = sgp4.FindPosition(current_time);
        CoordTopographic topo = obs.GetLookAngle(eci);

        /*
         * keep updating max elevation
         */
        if (topo.elevation > max_elevation)
        {
            max_elevation = topo.elevation;
        }
        else if (!fine_tune)
        {
            /*
             * passed max elevation
             * max elevation happened in the last 6 minutes
             * go back and fine tune max elevation value
             */
            current_time = current_time.AddSeconds(-2.0 * time_step);
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

        }
        else
        {
            /*
             * found max elevation
             */

            searching = false;
        }

        if (searching)
        {
            current_time = current_time.AddSeconds(time_step);
            if (current_time > los)
            {
                current_time = los;
            }
        }
    }

    return max_elevation;
}

DateTime FindCrossingPoint(
        const CoordGeodetic& user_geo,
        SGP4& sgp4,
        const DateTime& initial_time1,
        const DateTime& initial_time2,
        bool finding_aos)
{
    Observer obs(user_geo);

    bool searching = true;
    int cnt = 0;

    DateTime time1(initial_time1);
    DateTime time2(initial_time2);

    double diff = (time2 - time1).TotalSeconds();
    if (finding_aos)
    {
        diff = std::floor(diff);
    }
    else
    {
        diff = std::ceil(diff);
    }
    DateTime middle_time(time1.AddSeconds(diff));

    while (searching && cnt < 25)
    {
        /*
         * calculate satellite position
         */
        Eci eci = sgp4.FindPosition(middle_time);
        CoordTopographic topo = obs.GetLookAngle(eci);

        if (topo.elevation > 0.0)
        {
            if (finding_aos)
            {
                time2 = middle_time;
            }
            else
            {
                time1 = middle_time;
            }
        }
        else
        {
            if (finding_aos)
            {
                time1 = middle_time;
            }
            else
            {
                time2 = middle_time;
            }
        }

        /*
         * when two times are within a second, stop
         */
        if ((time2 - time1).TotalSeconds() < 1.5)
        {
            searching = false;
        }
        else
        {
            diff = (time2 - time1).TotalSeconds();
            if (finding_aos)
            {
                diff = std::floor(diff);
            }
            else
            {
                diff = std::ceil(diff);
            }
            middle_time = time1.AddSeconds(diff);
        }
        
        cnt++;
    }

    return middle_time;
}

std::list<struct PassDetails> GeneratePassList(
        const CoordGeodetic& user_geo,
        SGP4& sgp4,
        const DateTime& start_time,
        const DateTime& end_time,
        const int time_step)
{
    std::list<struct PassDetails> pass_list;

    Observer obs(user_geo);

    DateTime aos_time;
    DateTime los_time;

    bool found_aos = false;
    bool found_los = false;

    DateTime previous_time(start_time);
    DateTime current_time(start_time);

    while (current_time < end_time)
    {
        /*
         * calculate satellite position
         */
        Eci eci = sgp4.FindPosition(current_time);
        CoordTopographic topo = obs.GetLookAngle(eci);

        std::cout << std::fixed << current_time << " " << topo.elevation << std::endl;

        if (topo.elevation > 0.0)
        {
            /*
             * satellite is above horizon
             */
            if (start_time == current_time)
            {
                /*
                 * satellite was already above the horizon at the start time,
                 * so use the start time
                 */
                aos_time = start_time;
            }
            else
            {
                /*
                 * find the point at which the satellite crossed the horizon
                 */
                aos_time = FindCrossingPoint(
                        user_geo,
                        sgp4,
                        previous_time,
                        current_time,
                        true);
            }

            found_aos = true;
            found_los = false;
        }
        else if (found_aos)
        {
            /*
             * satellite now below horizon and we have an AOS, so record the
             * pass
             */
            los_time = FindCrossingPoint(
                    user_geo,
                    sgp4,
                    previous_time,
                    current_time,
                    false);

            found_aos = false;
            found_los = true;

            struct PassDetails pd;
            pd.aos = aos_time;
            pd.los = los_time;
            pd.max_elevation = FindMaxElevation(
                    user_geo,
                    sgp4,
                    aos_time,
                    los_time);

            pass_list.push_back(pd);
        }

        previous_time = current_time;

        if (found_los)
        {
            /*
             * at the end of the pass move the time along by 20mins
             */
            current_time = current_time + TimeSpan(0, 20, 0);
        }
        else
        {
            /*
             * move the time along by the time step value
             */
            current_time = current_time + TimeSpan(0, 0, time_step);
        }

        if (current_time > end_time)
        {
            /*
             * dont go past end time
             */
            current_time = end_time;
        }

        found_los = false;
    };

    if (found_aos)
    {
        /*
         * satellite still above horizon at end of search period, so use end
         * time as los
         */
        struct PassDetails pd;
        pd.aos = aos_time;
        pd.los = end_time;
        pd.max_elevation = FindMaxElevation(user_geo, sgp4, aos_time, los_time);
            
        pass_list.push_back(pd);
    }

    return pass_list;
}

int main() {

    CoordGeodetic geo(51.507406923983446, -0.12773752212524414, 0.05);
    Tle tle("UK-DMC 2                ",
        "1 25544U 98067A   12285.65009259  .00017228  00000-0  30018-3 0  4501",
        "2 25544 051.6477 262.7396 0017757 155.0745 185.1532 15.50683239796101");
    SGP4 sgp4(tle);

    std::cout << tle << std::endl;

    /*
     * generate 1 day schedule
     */
    DateTime start_date = DateTime::Now();
    DateTime end_date(start_date.AddDays(1.0));

    std::list<struct PassDetails> pass_list;
    std::list<struct PassDetails>::const_iterator itr;

    std::cout << "Start time: " << start_date << std::endl;
    std::cout << "End time  : " << end_date << std::endl << std::endl;

    /*
     * generate passes
     */
    pass_list = GeneratePassList(geo, sgp4, start_date, end_date, 180);

    if (pass_list.begin() == pass_list.end())
    {
        std::cout << "No passes found" << std::endl;
    }
    else
    {
        itr = pass_list.begin();
        do
        {
            std::cout
                << "AOS: " << itr->aos
                << ", LOS: " << itr->los
                << ", Max El: " << Util::RadiansToDegrees(itr->max_elevation)
                << std::endl;
        }
        while (++itr != pass_list.end());
    }

    return 0;
}
