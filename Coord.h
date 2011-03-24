#ifndef COORD_H_
#define COORD_H_

class CoordGeographic {
public:

    CoordGeographic()
    : lat_(0.0), lon_(0.0), alt_(0.0) {
    }

    CoordGeographic(double latitude, double longitude, double altitude)
    : lat_(latitude), lon_(longitude), alt_(altitude) {
    }

    virtual ~CoordGeographic() {
    };

    void SetLatitude(const double& latitude) {
        lat_ = latitude;
    }

    void SetLongitude(const double& longitude) {
        lon_ = longitude;
    }

    void SetAltitude(const double& altitude) {
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

class CoordTopographic {
public:

    CoordTopographic()
    : azimuth_(0.0), elevation_(0.0), range_(0.0), range_rate_(0.0) {
    }

    CoordTopographic(double azimuth, double elevation, double range, double range_rate)
    : azimuth_(azimuth), elevation_(elevation), range_(range), range_rate_(range_rate) {
    }

    virtual ~CoordTopographic() {
    };

    void SetAzimuth(const double& azimuth) {
        azimuth_ = azimuth;
    }

    void SetElevation(const double& elevation) {
        elevation_ = elevation;
    }

    void SetRange(const double& range) {
        range_ = range;
    }

    void SetRangeRate(const double& range_rate) {
        range_rate_ = range_rate;
    }

    double GetAzimuth() const {
        return azimuth_;
    }

    double GetElevation() const {
        return elevation_;
    }

    double GetRange() const {
        return range_;
    }

    double GetRangeRate() const {
        return range_rate_;
    }

private:
    /*
     * radians
     */
    double azimuth_;
    /*
     * radians
     */
    double elevation_;
    /*
     * kilometers
     */
    double range_;
    /*
     * kilometers / second
     */
    double range_rate_;
};

#endif

