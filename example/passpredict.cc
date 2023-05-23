#include "Topology.h"
#include<string.h>
int main()
{
    std::string path = "/mnt/c/Users/CProjects/handover/build/example";
    TopologyConstellation* topology = new TopologyConstellation(path);
    /*
    * generate 1 day schedule
    */
    libsgp4::DateTime start_date = libsgp4::DateTime::Now(true);
    libsgp4::DateTime end_date(start_date.AddHours(1.0));

    std::cout << "Start time: " << start_date << std::endl;
    std::cout << "End time  : " << end_date << std::endl << std::endl;

    /*
     * generate passes
     */
    topology->GeneratePassLists(start_date, end_date, 180);
    for (int i =0; i <topology->GetNumUserTerminals(); i++)
    {
        std::cout<<"-----------------------------------------city: "<<topology->GetUserEquipment((size_t)i)->GetName()<<" --------------------------------------------"<<std::endl;
        UserEquipment* ue = topology->GetUserEquipment((size_t)i);
        libsgp4::Observer obs = ue->GetObserver();
        libsgp4::DateTime time_now = start_date;
        while (time_now<end_date){
            Satellite * satellite_accessed = ue->GetSatelliteAccessed(time_now);
            if(satellite_accessed!=NULL) {
                std::cout << "time_now: " << time_now << " ue "<<ue->GetGid()<< " accesses satellite "
                          << ue->GetSatelliteAccessed(time_now)->GetSatelliteNumber()
                          << std::endl;
                /*
                * calculate satellite position
                 */
                libsgp4::Eci eci = satellite_accessed->GetSGP4Model().FindPosition(time_now);
                /*
                 * get look angle for observer to satellite
                 */
                libsgp4::CoordTopocentric topo = obs.GetLookAngle(eci);

                std::cout << topo << std::endl;
            }
            else{
                std::cout<<"No satellite can be accessed"<<std::endl;
            }
            time_now = time_now + libsgp4::TimeSpan(0,0,30);
        }
    }

    return 0;
}
