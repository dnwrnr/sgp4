#ifndef HANDOVER_TOPOLOGY_H
#define HANDOVER_TOPOLOGY_H
#include "Satellite.h"
#include "UserEquipment.h"
#include <vector>
#include <string>
#include <fstream>
#include<sstream>
#include <iostream>
class TopologyConstellation
{
public:

    TopologyConstellation(std::string constellation_network_dir);

    ~TopologyConstellation();

    int GetNumOrbits() const;
    int GetNumSatellitesPerOrbit() const;
    int GetNumSatellites() const;
    int GetNumUserTerminals() const;
    std::vector<Satellite *> GetSatellites() const;
    std::vector<UserEquipment *> GetUserEquipments() const;
    Satellite* GetSatellite(size_t sat_id) const;
    UserEquipment* GetUserEquipment(size_t ut_id) const;
    void GeneratePassLists(
            const libsgp4::DateTime& start_time,
            const libsgp4::DateTime& end_time,
            const int time_step);
private:
    // Build functions
    void ReadSatellites();
    void ReadUserTerminal();
    std::list<struct PassDetails>
    GeneratePassList(
            libsgp4::Observer& obs,
            libsgp4::SGP4& sgp4,
            const libsgp4::DateTime& start_time,
            const libsgp4::DateTime& end_time,
            const int time_step);
    double FindMaxElevation(
            libsgp4::Observer& obs,
            libsgp4::SGP4& sgp4,
            const libsgp4::DateTime& aos,
            const libsgp4::DateTime& los);
    libsgp4::DateTime FindCrossingPoint(
            libsgp4::Observer& obs,
            libsgp4::SGP4& sgp4,
            const libsgp4::DateTime& initial_time1,
            const libsgp4::DateTime& initial_time2,
            bool finding_aos);

    std::string m_constellation_network_dir;          //<! Directory containing satellite network information
    std::vector<Satellite *> m_satellites;           //<! Satellites
    std::vector<UserEquipment *> m_userEquipments;        //<! UserEquipments
    // Topology layout properties
    int m_num_satellites;
    int m_num_user_terminals;
    int m_num_orbits;
    int m_num_satellites_per_orbit;
};

#endif //HANDOVER_TOPOLOGY_H
