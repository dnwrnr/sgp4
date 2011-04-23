#ifndef COORDGEODETIC_H_
#define COORDGEODETIC_H_

class CoordGeodetic {
public:

    CoordGeodetic()
    : lat_(0.0), lon_(0.0), alt_(0.0) {
    }

    /*
     * radians
     */
    CoordGeodetic(double latitude, double longitude, double altitude)
    : lat_(latitude), lon_(longitude), alt_(altitude) {
    }

    CoordGeodetic(const CoordGeodetic& g);

    virtual ~CoordGeodetic() {
    };

    CoordGeodetic & operator =(const CoordGeodetic& b);
    bool operator ==(const CoordGeodetic& b) const;
    bool operator !=(const CoordGeodetic& b) const;

    void SetLatitude(const double latitude) {
        lat_ = latitude;
    }

    void SetLongitude(const double longitude) {
        lon_ = longitude;
    }

    void SetAltitude(const double altitude) {
        alt_ = altitude;
    }

    double GetLatitude() const {
        return lat_;
    }

    double GetLongitude() const {
        return lon_;
    }

    double GetAltitude() const {
        return alt_;
    }

private:
    /*
     * radians (north positive, south negative)
     */
    double lat_;
    /*
     * radians (east positive, west negative)
     */
    double lon_;
    /*
     * kilometers
     */
    double alt_;
};

#endif

