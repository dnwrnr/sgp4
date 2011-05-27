#ifndef COORDGEODETIC_H_
#define COORDGEODETIC_H_

#include <iostream>

struct CoordGeodetic {
public:

    CoordGeodetic()
    : latitude(0.0), longitude(0.0), altitude(0.0) {
    }

    /*
     * radians
     */
    CoordGeodetic(double lat, double lon, double alt)
    : latitude(lat), longitude(lon), altitude(alt) {
    }

    CoordGeodetic(const CoordGeodetic& g);

    virtual ~CoordGeodetic() {
    };

    CoordGeodetic & operator =(const CoordGeodetic& b);
    bool operator ==(const CoordGeodetic& b) const;
    bool operator !=(const CoordGeodetic& b) const;
    
    friend std::ostream& operator<< (std::ostream& stream, const CoordGeodetic& geo);
    
    /*
     * radians (north positive, south negative)
     */
    double latitude;
    /*
     * radians (east positive, west negative)
     */
    double longitude;
    /*
     * kilometers
     */
    double altitude;
};

#endif

