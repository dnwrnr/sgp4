#ifndef COORDTOPOGRAPHIC_H_
#define COORDTOPOGRAPHIC_H_

struct CoordTopographic {
public:

    CoordTopographic()
    : azimuth(0.0), elevation(0.0), range(0.0), range_rate(0.0) {
    }

    CoordTopographic(double az, double el, double rnge, double rnge_rate)
    : azimuth(az), elevation(el), range(rnge), range_rate(rnge_rate) {
    }

    CoordTopographic(const CoordTopographic& b);

    virtual ~CoordTopographic() {
    };

    CoordTopographic & operator =(const CoordTopographic& b);
    bool operator ==(const CoordTopographic& b) const;
    bool operator !=(const CoordTopographic& b) const;

    /*
     * radians
     */
    double azimuth;
    /*
     * radians
     */
    double elevation;
    /*
     * kilometers
     */
    double range;
    /*
     * kilometers / second
     */
    double range_rate;
};

#endif

