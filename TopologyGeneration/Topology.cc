#include "Topology.h"
#include "assert.h"
std::vector<std::string>
split_string(const std::string line, const std::string delimiter) {
    std::vector<std::string> result;
    std::string remainder = line;
    size_t idx = remainder.find(delimiter);
    while (idx != std::string::npos) {
        result.push_back(remainder.substr(0, idx));
        remainder = remainder.substr(idx + delimiter.size(), remainder.size());
        idx = remainder.find(delimiter);
    }
    result.push_back(remainder);
    return result;
}

/**
 * Split a string by the delimiter(s) and check that the split size is of expected size.
 * If it is not of expected size, throw an exception.
 *
 * @param line          Line (e.g., "a->b->c")
 * @param delimiter     Delimiter string (e.g., "->")
 * @param expected      Expected number (e.g., 3)
 *
 * @return Split vector (e.g., [a, b, c])
 */
std::vector<std::string> split_string_length(const std::string line, const std::string delimiter, size_t expected) {

    std::vector<std::string> the_split = split_string(line, delimiter);

    // It must match the expected split length, else throw exception
    if (the_split.size() != expected) {
        throw std::invalid_argument("invalid string split");
    }

    return the_split;
}

double TopologyConstellation::FindMaxElevation(
    libsgp4::Observer& obs,
    libsgp4::SGP4& sgp4,
    const libsgp4::DateTime& aos,
    const libsgp4::DateTime& los)
{
    bool running;

    double time_step = (los - aos).TotalSeconds() / 9.0;
    libsgp4::DateTime current_time(aos); //! current time
    libsgp4::DateTime time1(aos); //! start time of search period
    libsgp4::DateTime time2(los); //! end time of search period
    double max_elevation; //! max elevation

    running = true;

    do
    {
        running = true;
        max_elevation = -99999999999999.0;
        while (running && current_time < time2)
        {
            /*
             * find position
             */
            libsgp4::Eci eci = sgp4.FindPosition(current_time);
            libsgp4::CoordTopocentric topo = obs.GetLookAngle(eci);

            if (topo.elevation > max_elevation)
            {
                /*
                 * still going up
                 */
                max_elevation = topo.elevation;
                /*
                 * move time along
                 */
                current_time = current_time.AddSeconds(time_step);
                if (current_time > time2)
                {
                    /*
                     * dont go past end time
                     */
                    current_time = time2;
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
        time1 = current_time.AddSeconds(-2.0 * time_step);
        /*
         * make end time to current time
         */
        time2 = current_time;
        /*
         * current time to start time
         */
        current_time = time1;
        /*
         * recalculate time step
         */
        time_step = (time2 - time1).TotalSeconds() / 9.0;
    }
    while (time_step > 1.0);

    return max_elevation;
}

libsgp4::DateTime
TopologyConstellation::FindCrossingPoint(
    libsgp4::Observer& obs,
    libsgp4::SGP4& sgp4,
    const libsgp4::DateTime& initial_time1,
    const libsgp4::DateTime& initial_time2,
    bool finding_aos)
{
    bool running;
    int cnt;

    libsgp4::DateTime time1(initial_time1);
    libsgp4::DateTime time2(initial_time2);
    libsgp4::DateTime middle_time;

    running = true;
    cnt = 0;
    while (running && cnt++ < 16)
    {
        middle_time = time1.AddSeconds((time2 - time1).TotalSeconds() / 2.0);
        /*
         * calculate satellite position
         */
        libsgp4::Eci eci = sgp4.FindPosition(middle_time);
        libsgp4::CoordTopocentric topo = obs.GetLookAngle(eci);

        if (topo.elevation > 0.0)
        {
            /*
             * satellite above horizon
             */
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

        if ((time2 - time1).TotalSeconds() < 1.0)
        {
            /*
             * two times are within a second, stop
             */
            running = false;
            /*
             * remove microseconds
             */
            int us = middle_time.Microsecond();
            middle_time = middle_time.AddMicroseconds(-us);
            /*
             * step back into the pass by 1 second
             */
            middle_time = middle_time.AddSeconds(finding_aos ? 1 : -1);
        }
    }

    /*
     * go back/forward 1second until below the horizon
     */
    running = true;
    cnt = 0;
    while (running && cnt++ < 6)
    {
        libsgp4::Eci eci = sgp4.FindPosition(middle_time);
        libsgp4::CoordTopocentric topo = obs.GetLookAngle(eci);
        if (topo.elevation > 0)
        {
            middle_time = middle_time.AddSeconds(finding_aos ? -1 : 1);
        }
        else
        {
            running = false;
        }
    }

    return middle_time;
}

std::list<struct PassDetails>
TopologyConstellation::GeneratePassList(
        libsgp4::Observer& obs,
        libsgp4::SGP4& sgp4,
        const libsgp4::DateTime& start_time,
        const libsgp4::DateTime& end_time,
        const int time_step)
{
    std::list<struct PassDetails> pass_list;
    libsgp4::DateTime aos_time;
    libsgp4::DateTime los_time;

    bool found_aos = false;

    libsgp4::DateTime previous_time(start_time);
    libsgp4::DateTime current_time(start_time);

    while (current_time < end_time)
    {
        bool end_of_pass = false;

        /*
         * calculate satellite position
         */
        libsgp4::Eci eci = sgp4.FindPosition(current_time);
        libsgp4::CoordTopocentric topo = obs.GetLookAngle(eci);

        if (!found_aos && topo.elevation > 0.0)
        {
            /*
             * aos hasnt occured yet, but the satellite is now above horizon
             * this must have occured within the last time_step
             */
            if (start_time == current_time)
            {
                /*
                 * satellite was already above the horizon at the start,
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
                        obs,
                        sgp4,
                        previous_time,
                        current_time,
                        true);
            }
            found_aos = true;
        }
        else if (found_aos && topo.elevation < 0.0)
        {
            found_aos = false;
            /*
             * end of pass, so move along more than time_step
             */
            end_of_pass = true;
            /*
             * already have the aos, but now the satellite is below the horizon,
             * so find the los
             */
            los_time = FindCrossingPoint(
                    obs,
                    sgp4,
                    previous_time,
                    current_time,
                    false);

            struct PassDetails pd;
            pd.aos = aos_time;
            pd.los = los_time;
            pd.max_elevation = FindMaxElevation(
                    obs,
                    sgp4,
                    aos_time,
                    los_time);

            pass_list.push_back(pd);
        }

        /*
         * save current time
         */
        previous_time = current_time;

        if (end_of_pass)
        {
            /*
             * at the end of the pass move the time along by 30mins
             */
            current_time = current_time + libsgp4::TimeSpan(0, 30, 0);
        }
        else
        {
            /*
             * move the time along by the time step value
             */
            current_time = current_time + libsgp4::TimeSpan(0, 0, time_step);
        }

        if (current_time > end_time)
        {
            /*
             * dont go past end time
             */
            current_time = end_time;
        }
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
        pd.max_elevation = FindMaxElevation(obs, sgp4, aos_time, end_time);
        pass_list.push_back(pd);
    }

    return pass_list;
}

void
TopologyConstellation::GeneratePassLists(
    const libsgp4::DateTime& start_time,
    const libsgp4::DateTime& end_time,
    const int time_step)
{

    for (size_t i = 0; i < m_userEquipments.size(); ++i) {
        libsgp4::Observer obs(m_userEquipments.at(i)->GetObserver());
        std::list<struct PassDetails> pass_list_all;
        m_userEquipments.at(i)->SetStartAndEndTime(start_time,end_time);
        for (size_t j = 0; j < m_satellites.size(); ++j) {
            libsgp4::SGP4 sgp4(m_satellites.at(j)->GetSGP4Model());
            std::list<struct PassDetails> pass_list = GeneratePassList(obs,sgp4,start_time,end_time,time_step);
            if (pass_list.begin() != pass_list.end())
            {
                std::list<struct PassDetails>::const_iterator itr = pass_list.begin();
                do
                {
                    struct PassDetails new_detail;
                    new_detail.los = itr->los;
                    new_detail.aos = itr->aos;
                    new_detail.max_elevation = itr->max_elevation;
                    new_detail.satellite = m_satellites.at(j);
                    pass_list_all.push_back(new_detail);
                }while (++itr != pass_list.end());
            }
        }
        m_userEquipments.at(i)->SetPassList(pass_list_all);
    }
    for (size_t i = 0; i < m_userEquipments.size(); ++i) {
        m_userEquipments.at(i)->GenerateAccessList();
    }
}

TopologyConstellation::TopologyConstellation(std::string constellation_network_dir)
{
    m_constellation_network_dir = constellation_network_dir;
    ReadSatellites();
    ReadUserTerminal();
}

TopologyConstellation::~TopologyConstellation()
{
    //nothing
}

void TopologyConstellation::ReadSatellites()
{
    std::cout << "  > Read TLE file of satellites" << std::endl;
    // Open file
    std::ifstream fs;
    fs.open(m_constellation_network_dir + "/tles.txt");
    assert(fs.is_open());

    // First line:
    // <orbits> <satellites per orbit>
    std::string orbits_and_n_sats_per_orbit;
    std::getline(fs, orbits_and_n_sats_per_orbit);
    std::vector<std::string> res = split_string_length(orbits_and_n_sats_per_orbit, " ", 2);
    m_num_orbits = std::stoi(res[0], 0, 10);
    m_num_satellites_per_orbit = std::stoi(res[1], 0, 10);
    int num_satellites = m_num_orbits * m_num_satellites_per_orbit;

    // Associate satellite mobility model with each node
    uint32_t counter = 0;
    std::string name, tle1, tle2;
    while (std::getline(fs, name)) {
        std::getline(fs, tle1);
        std::getline(fs, tle2);

        // Format:
        // <name>
        // <TLE line 1>
        // <TLE line 2>

        // Create satellite
        Satellite* satellite = new Satellite();
        satellite->SetName(name);
        satellite->SetTleInfo(tle1, tle2);
        satellite->SetSatelliteNumber((int)counter);

        // Add to all satellites present
        m_satellites.push_back(satellite);

        counter++;
    }

    // Check that exactly that number of satellites has been read in
    if ((int)counter != num_satellites) {
        throw std::runtime_error("Number of satellites defined in the TLEs does not match");
    }


    fs.close();
    std::cout << "  > Number of satellites........ " << num_satellites << std::endl;
    m_num_satellites = num_satellites;
}

void TopologyConstellation::ReadUserTerminal()
{
    std::cout << "  > Read user terminals" << std::endl;
    // Create a new file stream to open the file
    std::ifstream fs;
    fs.open(m_constellation_network_dir + "/user_terminals.txt");
    assert(fs.is_open());
    int count_ue = 0;
    // Read user terminals from each line
    std::string line;
    while (std::getline(fs, line)) {
        std::vector<std::string> res = split_string_length(line, ",", 5);
        // All eight values
        int gid = std::stoi(res[0], 0, 10);
        std::string name = res[1];
        double latitude = std::stod(res[2]);
        double longitude = std::stod(res[3]);
        double elevation = std::stod(res[4]);
        UserEquipment * ue = new UserEquipment((uint32_t)gid, name, latitude, longitude, elevation);
        m_userEquipments.push_back(ue);
        count_ue ++;
    }

    fs.close();
    std::cout << "  > Number of user terminals... " << m_userEquipments.size() << std::endl;
    assert(count_ue == (int)(m_userEquipments.size()));
    m_num_user_terminals = count_ue;
}

int
TopologyConstellation::GetNumOrbits() const
{
    return m_num_orbits;
}
int
TopologyConstellation::GetNumSatellitesPerOrbit() const
{
    return m_num_satellites_per_orbit;
}
int
TopologyConstellation::GetNumSatellites() const
{
    return m_num_satellites;
}
int
TopologyConstellation::GetNumUserTerminals() const
{
    return m_num_user_terminals;
}

std::vector<Satellite *>
TopologyConstellation::GetSatellites() const
{
    return m_satellites;
}

Satellite *
TopologyConstellation::GetSatellite(size_t sat_id) const
{
    assert(sat_id<(size_t)m_num_satellites);
    return m_satellites.at(sat_id);
}


std::vector<UserEquipment *>
TopologyConstellation::GetUserEquipments() const
{
    return m_userEquipments;
}

UserEquipment*
TopologyConstellation::GetUserEquipment(size_t ut_id) const
{
    assert(ut_id<(size_t)m_num_user_terminals);
    return m_userEquipments.at(ut_id);
}
