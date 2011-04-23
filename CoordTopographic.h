#ifndef COORDTOPOGRAPHIC_H_
#define COORDTOPOGRAPHIC_H_

class CoordTopographic {
public:

    CoordTopographic()
    : azimuth_(0.0), elevation_(0.0), range_(0.0), range_rate_(0.0) {
    }

    CoordTopographic(double azimuth, double elevation, double range, double range_rate)
    : azimuth_(azimuth), elevation_(elevation), range_(range), range_rate_(range_rate) {
    }

    CoordTopographic(const CoordTopographic& b);

    virtual ~CoordTopographic() {
    };

    CoordTopographic & operator =(const CoordTopographic& b);
    bool operator ==(const CoordTopographic& b) const;
    bool operator !=(const CoordTopographic& b) const;

    void SetAzimuth(const double azimuth) {
        azimuth_ = azimuth;
    }

    void SetElevation(const double elevation) {
        elevation_ = elevation;
    }

    void SetRange(const double range) {
        range_ = range;
    }

    void SetRangeRate(const double range_rate) {
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

