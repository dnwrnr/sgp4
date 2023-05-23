#include "Satellite.h"
#include "assert.h"

Satellite::Satellite (void) :
        m_name (""), m_tle1 (""), m_tle2 ("")
{

}

bool
Satellite::IsInitialized (void) const
{
    return ((m_tle1 != "") && (m_tle2 != ""));
}

void
Satellite::SetName (const std::string &name)
{
    m_name = name.substr (0, name.find_last_not_of (" ") + 1);
}

int
Satellite::GetSatelliteNumber (void) const
{
    return (IsInitialized () ? m_satellite_Id : 0);
}

void
Satellite::SetSatelliteNumber(int satellite_Id)
{
    m_satellite_Id = satellite_Id;
}

std::string
Satellite::GetName (void) const
{
    return m_name;
}

std::pair<std::string, std::string>
Satellite::GetTleInfo (void) const
{
    return std::make_pair (m_tle1, m_tle2);
}


libsgp4::Eci
Satellite::GetEci(const libsgp4::DateTime &dt) const
{
    assert(IsInitialized ());
    return m_satellite->FindPosition(dt);
}

libsgp4::SGP4
Satellite::GetSGP4Model() const
{
    assert(IsInitialized ());
    return (*m_satellite);
}

libsgp4::CoordGeodetic
Satellite::GetCoordGeodetic(const libsgp4::DateTime& dt) const
{
    if (!IsInitialized ())
        return libsgp4::CoordGeodetic();
    return (m_satellite->FindPosition(dt)).ToGeodetic();
}

libsgp4::DateTime
Satellite::GetTleEpoch (void) const {
    if (!IsInitialized ())
        return libsgp4::DateTime ();

    return m_epoch;
}

bool
Satellite::SetTleInfo (const std::string &line1, const std::string &line2)
{
    m_tle1 = std::string (line1.c_str ());
    m_tle2 = std::string (line2.c_str ());
    size_t pos = m_tle1.find('\r');
    if (pos != std::string::npos) {
        m_tle1.erase(pos, 1);
    }
    pos = m_tle2.find('\r');
    if (pos != std::string::npos) {
        m_tle2.erase(pos, 1);
    }
    libsgp4::Tle tle(m_name,
                     m_tle1,
                     m_tle2);
    m_satellite = new libsgp4::SGP4(tle);
    m_epoch = tle.Epoch();
    return true;
}






