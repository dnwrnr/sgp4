/*
 * Copyright 2013 Daniel Warner <contact@danrw.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include <libsgp4/CoordGeodetic.h>
#include <libsgp4/CoordTopocentric.h>
#include <libsgp4/Observer.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Util.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>

struct PassDetails
{
    libsgp4::DateTime aos;
    libsgp4::DateTime los;
    double maxElevation;
};

double FindMaxElevation(
        const libsgp4::CoordGeodetic& userGeo,
        libsgp4::SGP4& sgp4,
        const libsgp4::DateTime& aos,
        const libsgp4::DateTime& los)
{
    libsgp4::Observer obs(userGeo);

    bool running;

    double timeStep = (los - aos).TotalSeconds() / 9.0;
    libsgp4::DateTime currentTime(aos); //! current time
    libsgp4::DateTime time1(aos); //! start time of search period
    libsgp4::DateTime time2(los); //! end time of search period
    double maxElevation; //! max elevation

    running = true;

    do
    {
        running = true;
        maxElevation = -99999999999999.0;
        while (running && currentTime < time2)
        {
            /*
             * find position
             */
            libsgp4::Eci eci = sgp4.FindPosition(currentTime);
            libsgp4::CoordTopocentric topo = obs.GetLookAngle(eci);

            if (topo.elevation > maxElevation)
            {
                /*
                 * still going up
                 */
                maxElevation = topo.elevation;
                /*
                 * move time along
                 */
                currentTime = currentTime.AddSeconds(timeStep);
                if (currentTime > time2)
                {
                    /*
                     * dont go past end time
                     */
                    currentTime = time2;
                }
            }
            else
            {
                /*
                 * stop
                 */
                running = false;
            }
        }

        /*
         * make start time to 2 time steps back
         */
        time1 = currentTime.AddSeconds(-2.0 * timeStep);
        /*
         * make end time to current time
         */
        time2 = currentTime;
        /*
         * current time to start time
         */
        currentTime = time1;
        /*
         * recalculate time step
         */
        timeStep = (time2 - time1).TotalSeconds() / 9.0;
    }
    while (timeStep > 1.0);

    return maxElevation;
}

libsgp4::DateTime FindCrossingPoint(
        const libsgp4::CoordGeodetic& userGeo,
        libsgp4::SGP4& sgp4,
        const libsgp4::DateTime& initialTime1,
        const libsgp4::DateTime& initialTime2,
        bool findingAos)
{
    libsgp4::Observer obs(userGeo);

    bool running;
    int cnt;

    libsgp4::DateTime time1(initialTime1);
    libsgp4::DateTime time2(initialTime2);
    libsgp4::DateTime middleTime;

    running = true;
    cnt = 0;
    while (running && cnt++ < 16)
    {
        middleTime = time1.AddSeconds((time2 - time1).TotalSeconds() / 2.0);
        /*
         * calculate satellite position
         */
        libsgp4::Eci eci = sgp4.FindPosition(middleTime);
        libsgp4::CoordTopocentric topo = obs.GetLookAngle(eci);

        if (topo.elevation > 0.0)
        {
            /*
             * satellite above horizon
             */
            if (findingAos)
            {
                time2 = middleTime;
            }
            else
            {
                time1 = middleTime;
            }
        }
        else
        {
            if (findingAos)
            {
                time1 = middleTime;
            }
            else
            {
                time2 = middleTime;
            }
        }

        if ((time2 - time1).TotalSeconds() < 1.0)
        {
            /*
             * two times are within a second, stop
             */
            running = false;
            /*
             * remove microseconds
             */
            int us = middleTime.Microsecond();
            middleTime = middleTime.AddMicroseconds(-us);
            /*
             * step back into the pass by 1 second
             */
            middleTime = middleTime.AddSeconds(findingAos ? 1 : -1);
        }
    }

    /*
     * go back/forward 1second until below the horizon
     */
    running = true;
    cnt = 0;
    while (running && cnt++ < 6)
    {
        libsgp4::Eci eci = sgp4.FindPosition(middleTime);
        libsgp4::CoordTopocentric topo = obs.GetLookAngle(eci);
        if (topo.elevation > 0)
        {
            middleTime = middleTime.AddSeconds(findingAos ? -1 : 1);
        }
        else
        {
            running = false;
        }
    }

    return middleTime;
}

std::list<struct PassDetails> GeneratePassList(
        const libsgp4::CoordGeodetic& userGeo,
        libsgp4::SGP4& sgp4,
        const libsgp4::DateTime& startTime,
        const libsgp4::DateTime& endTime,
        const int timeStep)
{
    std::list<struct PassDetails> passList;

    libsgp4::Observer obs(userGeo);

    libsgp4::DateTime aosTime;
    libsgp4::DateTime losTime;

    bool foundAos = false;

    libsgp4::DateTime previousTime(startTime);
    libsgp4::DateTime currentTime(startTime);

    while (currentTime < endTime)
    {
        bool endOfPass = false;

        /*
         * calculate satellite position
         */
        libsgp4::Eci eci = sgp4.FindPosition(currentTime);
        libsgp4::CoordTopocentric topo = obs.GetLookAngle(eci);

        if (!foundAos && topo.elevation > 0.0)
        {
            /*
             * aos hasnt occured yet, but the satellite is now above horizon
             * this must have occured within the last timeStep
             */
            if (startTime == currentTime)
            {
                /*
                 * satellite was already above the horizon at the start,
                 * so use the start time
                 */
                aosTime = startTime;
            }
            else
            {
                /*
                 * find the point at which the satellite crossed the horizon
                 */
                aosTime = FindCrossingPoint(
                        userGeo,
                        sgp4,
                        previousTime,
                        currentTime,
                        true);
            }
            foundAos = true;
        }
        else if (foundAos && topo.elevation < 0.0)
        {
            foundAos = false;
            /*
             * end of pass, so move along more than timeStep
             */
            endOfPass = true;
            /*
             * already have the aos, but now the satellite is below the horizon,
             * so find the los
             */
            losTime = FindCrossingPoint(
                    userGeo,
                    sgp4,
                    previousTime,
                    currentTime,
                    false);

            struct PassDetails pd;
            pd.aos = aosTime;
            pd.los = losTime;
            pd.maxElevation = FindMaxElevation(
                    userGeo,
                    sgp4,
                    aosTime,
                    losTime);

            passList.push_back(pd);
        }

        /*
         * save current time
         */
        previousTime = currentTime;

        if (endOfPass)
        {
            /*
             * at the end of the pass move the time along by 30mins
             */
            currentTime = currentTime + libsgp4::TimeSpan(0, 30, 0);
        }
        else
        {
            /*
             * move the time along by the time step value
             */
            currentTime = currentTime + libsgp4::TimeSpan(0, 0, timeStep);
        }

        if (currentTime > endTime)
        {
            /*
             * dont go past end time
             */
            currentTime = endTime;
        }
    };

    if (foundAos)
    {
        /*
         * satellite still above horizon at end of search period, so use end
         * time as los
         */
        struct PassDetails pd;
        pd.aos = aosTime;
        pd.los = endTime;
        pd.maxElevation = FindMaxElevation(userGeo, sgp4, aosTime, endTime);
        passList.push_back(pd);
    }

    return passList;
}

int main()
{
    libsgp4::CoordGeodetic geo(51.507406923983446, -0.12773752212524414, 0.05);
    libsgp4::Tle tle("GALILEO-PFM (GSAT0101)  ",
        "1 37846U 11060A   12293.53312491  .00000049  00000-0  00000-0 0  1435",
        "2 37846  54.7963 119.5777 0000994 319.0618  40.9779  1.70474628  6204");
    libsgp4::SGP4 sgp4(tle);

    std::cout << tle << std::endl;

    libsgp4::DateTime startDate = libsgp4::DateTime::Now(true);
    libsgp4::DateTime endDate(startDate.AddDays(7.0));

    std::cout << "Start time: " << startDate << std::endl;
    std::cout << "End time  : " << endDate << std::endl << std::endl;

    try
    {
        std::list<struct PassDetails> passList = GeneratePassList(geo, sgp4, startDate, endDate, 180);

        if (passList.begin() == passList.end())
        {
            std::cout << "No passes found" << std::endl;
        }
        else
        {
            std::stringstream ss;

            ss << std::right << std::setprecision(1) << std::fixed;

            std::list<struct PassDetails>::const_iterator itr = passList.begin();
            do
            {
                ss  << "AOS: " << itr->aos
                    << ", LOS: " << itr->los
                    << ", Max El: " << std::setw(4) << libsgp4::Util::RadiansToDegrees(itr->maxElevation)
                    << ", Duration: " << (itr->los - itr->aos)
                    << std::endl;
            }
            while (++itr != passList.end());

            std::cout << ss.str();
        }
    }
    catch (libsgp4::SatelliteException& e)
    {
        std::cerr << "Satellite error: " << e.what() << std::endl;
        return 1;
    }
    catch (libsgp4::DecayedException& e)
    {
        std::cerr << "Satellite decayed: " << e.what() << std::endl;
        std::cerr << "Position: " << e.Position() << std::endl;
        return 1;
    }

    return 0;
}
