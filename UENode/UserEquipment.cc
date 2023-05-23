#include "UserEquipment.h"
#include "assert.h"
#include "Util.h"
UserEquipment::UserEquipment(uint32_t gid, std::string name, double latitude, double longitude, double altitude)
{
    m_gid = gid;
    m_name = name;
    m_latitude = latitude;
    m_longitude = longitude;
    m_altitude =altitude;
    libsgp4::CoordGeodetic user_geo(m_latitude, m_longitude, m_altitude);
    m_observer = new libsgp4::Observer(user_geo);
    m_elevation_threshold = 10.0;
}

UserEquipment::~UserEquipment()
{

}


uint32_t
UserEquipment::GetGid()
{
    return m_gid;
}

std::string
UserEquipment::GetName()
{
    return m_name;
}

double
UserEquipment::GetLatitude()
{
    return m_latitude;
}

double
UserEquipment::GetLongitude()
{
    return m_longitude;
}

double
UserEquipment::GetAltitude()
{
    return m_altitude;
}

void
UserEquipment::SetPassList(PassList passList)
{
    m_passList = passList;
}

PassList
UserEquipment::GetPassList()
{
    return m_passList;
}



Satellite*
UserEquipment::GetSatelliteAccessed(libsgp4::DateTime dt){
    std::list<AccessDetails>::const_iterator iter = m_accessList.begin();
    Satellite * satellite_accessed = NULL;
    while(iter!=m_accessList.end()){
        if (dt >= iter->aos && dt<=iter->los)
        {
            satellite_accessed = iter->satellite;
            break;
        }
        iter++;
    }
    return satellite_accessed;
}


libsgp4::Observer
UserEquipment::GetObserver() const
{
    return (*m_observer);
}

void
UserEquipment::SetStartAndEndTime(libsgp4::DateTime start_time, libsgp4::DateTime end_time)
{
    m_start_time = start_time;
    m_end_time = end_time;
}

struct AccessDetails
UserEquipment::AccessSatellite(libsgp4::DateTime dt)
{
    std::list<struct PassDetails> passList = GetPassList();
    libsgp4::DateTime time_now = dt;
    double max_elevation = 0.0;
    Satellite* satellite_accessed = NULL;
    libsgp4::DateTime access_time;
    libsgp4::DateTime lose_time;
    struct AccessDetails ad;
    std::list<struct PassDetails>::const_iterator itr = passList.begin();
    do
    {
        access_time = itr->aos;
        lose_time = itr->los;
        double elevation = 0.0;
        if(time_now > lose_time || time_now < access_time){
            continue;
        }else{
            /*
            * calculate satellite position
            */
            libsgp4::Eci eci = itr->satellite->GetSGP4Model().FindPosition(dt);
            /*
             * get look angle for observer to satellite
             */
            elevation = libsgp4::Util::RadiansToDegrees(m_observer->GetLookAngle(eci).elevation);
            if(max_elevation < elevation)
            {
                max_elevation = elevation;
                access_time = itr->aos;
                lose_time = itr->los;
                satellite_accessed = itr->satellite;
            }
        }
    }while (++itr != passList.end());
    ad.satellite = satellite_accessed;
    ad.aos = time_now;
    ad.los = lose_time;
    return ad;
}

void
UserEquipment::GenerateAccessList()
{
    libsgp4::DateTime time_now = m_start_time;
    Satellite* satellite_accessed;
    AccessDetails ad;
    //initial access
    while (time_now < m_end_time) {
        ad = AccessSatellite(time_now);
        if (ad.satellite != NULL) {
            satellite_accessed = ad.satellite;
            break;
        }
        time_now = time_now + libsgp4::TimeSpan(0, 0, 10);
    }
    //handover
    while (time_now < m_end_time){
        time_now = time_now + libsgp4::TimeSpan(0, 0, 10);
        libsgp4::Eci eci = satellite_accessed->GetSGP4Model().FindPosition(time_now);
        /*
         * get look angle for observer to satellite
         */
        double elevation = libsgp4::Util::RadiansToDegrees(m_observer->GetLookAngle(eci).elevation);
        if(elevation < m_elevation_threshold||time_now>=ad.los)
        {
            //handover
            ad.los = time_now;
            m_accessList.push_back(ad);
            while (time_now < m_end_time) {
                time_now = time_now + libsgp4::TimeSpan(0, 0, 10);
                ad = AccessSatellite(time_now);
                if (ad.satellite != NULL) {
                    satellite_accessed = ad.satellite;
                    break;
                }
            }
        }
    }
}

