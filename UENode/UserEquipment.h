#ifndef HANDOVER_USEREQUIPMENT_H
#define HANDOVER_USEREQUIPMENT_H
#include <Observer.h>
#include <SGP4.h>
#include <Util.h>
#include <CoordTopocentric.h>
#include <CoordGeodetic.h>
#include <Satellite.h>
#include <cmath>
#include <iostream>
#include <list>
#include "Satellite.h"
class Satellite;
struct PassDetails
{
    libsgp4::DateTime aos;
    libsgp4::DateTime los;
    double max_elevation;
    Satellite* satellite;
};

struct AccessDetails
{
    libsgp4::DateTime aos;
    libsgp4::DateTime los;
    Satellite* satellite;
};

typedef std::list<PassDetails> PassList;
typedef std::list<AccessDetails> AccessList;

class Satellite;
class UserEquipment
{
public:

    /**
    * @brief Constructor function.
    * @param gid the UE's number (ID).
    * @param name the UE's name.
    * @param latitude the UE's latitude in degree °.
    * @param longitude the UE's longitude in degree °.
    * @param altitude the UE's altitude in km.
    */
    UserEquipment(
            uint32_t gid, std::string name,
            double latitude, double longitude, double altitude
    );
    ~UserEquipment();

    /**
    * @brief Retrieve the UE's number (UE_ID).
    * @return the UE's number (ID).
    */
    uint32_t GetGid();

    /**
     * @brief Retrieve UE's name.
     * @return the UE's name.
     */
    std::string GetName();

    /**
     * @brief Retrieve UE's Latitude.
     * @return the UE's Latitude.
     */
    double GetLatitude();

    /**
     * @brief Retrieve UE's Longitude.
     * @return the UE's Longitude.
     */
    double GetLongitude();

    /**
     * @brief Retrieve UE's Longitude.
     * @return the UE's Longitude.
     */
    double GetAltitude();

    /**
     * @brief Set UE's ephemeris.
     * @param passList ephemeris.
     */
    void SetPassList(PassList passList);

    /**
     * @brief Retrieve UE's ephemeris.
     * @return the UE's ephemeris.
     */
    PassList GetPassList();

    /**
     * @brief Get current UE's satellite Id has been accessed.
     * @return the satellite Id.
     */
    int GetSatelliteIdAccessed();

    /**
     * @brief Set satellite Id of UE will be accessed.
     * @param satelliteId the satellite Id.
     */
    void SetSatelliteIdAccessed(int satelliteId);

    /**
     * @brief Get current UE's ptr<satellite> has been accessed.
     * @param dt time
     * @return the satellite.
     */
    Satellite* GetSatelliteAccessed(libsgp4::DateTime dt);


    /**
     * @brief Get current UE's libsgp4::Observer.
     * @return the observer.
     */
    libsgp4::Observer GetObserver() const;

    /**
     * @brief Access the satellite with maximum inclination
     * @param dt delay from simulation start.
     * @return The satellite with maximum elevation will be accessed.
     *
     */
    struct AccessDetails AccessSatellite (libsgp4::DateTime dt);

    /**
     * @brief Generate access list
     */
    void GenerateAccessList();

    /**
     * @brief Set start time and end time of simulation.
     *
     */
    void SetStartAndEndTime(libsgp4::DateTime start_time, libsgp4::DateTime end_time);

private:
    uint32_t m_gid;        // Unique UTs identifier
    std::string m_name;    // Name
    double m_latitude;     // Latitude degree °
    double m_longitude;    // Longitude degree °
    double m_altitude;    // Altitude meter m
    PassList m_passList; // passList
    AccessList m_accessList; //accessList
    libsgp4::Observer *m_observer;
    libsgp4::DateTime m_start_time;
    libsgp4::DateTime m_end_time;
    double m_elevation_threshold;

};

#endif //HANDOVER_USEREQUIPMENT_H
