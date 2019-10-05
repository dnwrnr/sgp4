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

#include <CoordTopocentric.h>
#include <CoordGeodetic.h>
#include <Observer.h>
#include <SGP4.h>

#include <iostream>
#include <getopt.h>
#include <curl/curl.h>

void usage()
{
    printf("isstrack [options]\n");
    printf("  h     : print usage\n");
    printf("  t   : Observer Latitude\n");
    printf("  g   : Observer Longitude\n");
    printf("  a   : Observer Altitude (km)\n");
}

size_t write_data(char *ptr, size_t size, size_t nmemb, void *userdata) {
    std::ostringstream *stream = (std::ostringstream*)userdata;
    size_t count = size * nmemb;
    stream->write(ptr, count);
    return count;
}

int main(int argc, char*argv[])
{
    // varialbles
    std::stringstream stream;  //stringstream with TLEs from Celestrack
    std::string name;
    std::string line1;
    std::string line2;

    // options
    float lat = 51.507406923983446;  // Latitude
    float lng = -0.127737522125244;  // Longitude
    float alt = 0.05;                // Altitude (km)

    int opt;
    while((opt = getopt(argc,argv,"hn:t:g:a:")) != EOF){
        switch (opt) {
        case 'h': usage();                      return 0;
        case 't': lat = atof(optarg);   break;
        case 'g': lng = atof(optarg);   break;
        case 'a': alt = atof(optarg);   break;
        default:
            exit(1);
        }
    }

    // configure and execute lib curl
    CURL *curl;
    curl = curl_easy_init();
    curl_easy_setopt(curl, CURLOPT_URL, "http://www.celestrak.com/norad/elements/stations.txt");
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &stream);
    curl_easy_perform(curl);
    curl_easy_cleanup(curl);

    // get ISS TLE from stringstream
    std::string line;
    while(std::getline(stream,line)){
      if(line.find("ISS (ZARYA)") != std::string::npos)
      {
        // printf("%s\n",line.c_str());
        name = line.c_str();
        getline(stream,line);
        line.pop_back();
        line1 = line.c_str();
        // printf("%s\n",line.c_str());
        getline(stream,line);
        line.pop_back();
        line2 = line.c_str();
        // printf("%s\n",line.c_str());
      }
    }

    // SGP4 objects
    Tle tle = Tle(name, line1, line2);
    SGP4 sgp4(tle);
    Observer obs(lat, lng, alt);

    DateTime ft = DateTime::Now(true);    //future time
    while (true) {
      DateTime dt = DateTime::Now(true);

      /*
       * calculate satellite position
       */
      Eci eci = sgp4.FindPosition(dt);
      /*
       * convert satellite position to geodetic coordinates
       */
      CoordGeodetic geo = eci.ToGeodetic();

      if (dt > ft) {
        ft = DateTime::Now(true).AddSeconds(1);
        std::cout << dt << " Current Position: " << geo << '\n';
      }
    }

    return 0;
}
