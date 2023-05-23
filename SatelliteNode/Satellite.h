#ifndef HANDOVER_SATELLITE_H
#define HANDOVER_SATELLITE_H
#include <SGP4.h>
#include <Util.h>
#include <CoordTopocentric.h>
#include <CoordGeodetic.h>
#include <Vector.h>
#include <cmath>
#include <iostream>
#include <list>
#include <UserEquipment.h>

class UserEquipment;
class Satellite
{
public:
    Satellite();
    ~Satellite();

    /**
   * @brief Retrieve the satellite's number (SAT_ID).
   * @return the satellite's number (ID) if it has already been
   *         initialized or 0 otherwise.
   */
    int GetSatelliteNumber (void) const;

    /**
    * @brief Set the satellite's number (NORAD SAT_ID).
    * @param the satellite's unique Id.
    */
    void SetSatelliteNumber(int satellite_Id);

    /**
     * @brief Retrieve satellite's name.
     * @return the satellite's name or an empty string if has not yet been set.
     */
    std::string GetName (void) const;

    /**
     * @brief Retrieve the TLE information used to initialize this satellite.
     * @return an std::pair with the two lines used to initialize this satellite,
     *         or an std::pair with two empty strings if it has not yet been set.
     */
    std::pair<std::string, std::string> GetTleInfo (void) const;


    /**
     * @brief Retrieve the TLE epoch time.
     * @return the TLE epoch time or 0001/01/01 00:00:00.000000 if the satellite has not
     *         yet been initialized.
     */
    libsgp4::DateTime GetTleEpoch (void) const;

    /**
   * @brief Set satellite's name.
   * @param name Satellite's name.
   */
    void SetName (const std::string &name);

    /**
     * @brief Set satellite's TLE information required for its initialization.
     * @param line1 First line of the TLE data format.
     * @param line2 Second line of the TLE data format.
     * @param name Satellite's name.
     * @return a boolean indicating whether the initialization succeeded.
     */
    bool SetTleInfo (const std::string &line1, const std::string &line2);

    /**
   * @brief Get the prediction for the satellite's position at a given time.
   * @param dt When.
   * @return a libsgp4::Eci object containing the satellite's position
   */
    libsgp4::Eci GetEci(const libsgp4::DateTime& dt) const;


    /**
     * @brief Get the prediction for the satellite's Earth-centered inertial position at a given time.
     * @param dt When.
     * @return CoordGeodetic of satellite latitude, longitude, altitude
     */
    libsgp4::CoordGeodetic GetCoordGeodetic(const libsgp4::DateTime& dt) const;


    /**
     * @brief Get current satellie's libsgp4::SGP4.
     * @return the libsgp4::SGP4.
     */
    libsgp4::SGP4 GetSGP4Model() const;

private:
    /**
     * @brief Check if the satellite has already been initialized.
     * @return a boolean indicating whether the satellite is initialized.
     */
    bool IsInitialized (void) const;

private:
    std::string m_name;                               //!< satellite's name.
    std::string m_tle1, m_tle2;                       //!< satellite's TLE data.
    libsgp4::SGP4 * m_satellite;
    int m_satellite_Id;
    libsgp4::DateTime m_epoch;
};
#endif //HANDOVER_SATELLITE_H
