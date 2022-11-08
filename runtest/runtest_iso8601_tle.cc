/*
 * Copyright 2022 Andy Kirkham
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


#include <Tle.h>
#include <DateTime.h>

#include <string>

/*
Example of current TLE format in JSON from Celestrak.
[{
    "OBJECT_NAME": "STARLINK-1007",
    "OBJECT_ID": "2019-074A",
    "EPOCH": "2022-11-08T06:14:56.037120",
    "MEAN_MOTION": 15.06405436,
    "ECCENTRICITY": 0.0001911,
    "INCLINATION": 53.0559,
    "RA_OF_ASC_NODE": 251.8795,
    "ARG_OF_PERICENTER": 48.9031,
    "MEAN_ANOMALY": 311.2123,
    "EPHEMERIS_TYPE": 0,
    "CLASSIFICATION_TYPE": "U",
    "NORAD_CAT_ID": 44713,
    "ELEMENT_SET_NO": 999,
    "REV_AT_EPOCH": 16525,
    "BSTAR": 0.00033293,
    "MEAN_MOTION_DOT": 4.682e-5,
    "MEAN_MOTION_DDOT": 0
}]
*/

int RunTest(void)
{
    TleArgs args;
    args.name = std::string("STARLINK-1007");
    args.int_designator = std::string("2019-074A");
    args.epoch = std::string("2022-11-08T06:14:56.037120");
    args.ephemeris_type = 0;
    args.classification_type = std::string("U");
    args.mean_motion = 15.06405436;
    args.mean_anomaly = 311.2123;
    args.inclination = 53.0559;
    args.right_ascending_node = 251.8795;
    args.eccentricity = 0.0001911;
    args.argument_perigee = 48.9031;
    args.bstar = 0.00033293;
    args.norad_number = 44713;
    args.orbit_number = 16525;
    args.mean_motion_dot = 4.682e-5;
    args.mean_motion_ddot = 0;

    Tle dt(args);
    std::cout << dt.ToString();
    return 1;
}

int main()
{
    return RunTest();
}
