#ifndef ORBITALELEMENTS_H_
#define ORBITALELEMENTS_H_

#include "Tle.h"

class OrbitalElements
{
public:
    OrbitalElements(const Tle& tle);

    virtual ~OrbitalElements()
    {
    }

    /*
     * XMO
     */
    double MeanAnomoly() const
    {
        return mean_anomoly_;
    }

    /*
     * XNODEO
     */
    double AscendingNode() const
    {
        return ascending_node_;
    }

    /*
     * OMEGAO
     */
    double ArgumentPerigee() const
    {
        return argument_perigee_;
    }

    /*
     * EO
     */
    double Eccentricity() const
    {
        return eccentricity_;
    }

    /*
     * XINCL
     */
    double Inclination() const
    {
        return inclination_;
    }

    /*
     * XNO
     */
    double MeanMotion() const
    {
        return mean_motion_;
    }

    /*
     * BSTAR
     */
    double BStar() const
    {
        return bstar_;
    }

    /*
     * AODP
     */
    double RecoveredSemiMajorAxis() const
    {
        return recovered_semi_major_axis_;
    }

    /*
     * XNODP
     */
    double RecoveredMeanMotion() const
    {
        return recovered_mean_motion_;
    }

    /*
     * PERIGE
     */
    double Perigee() const
    {
        return perigee_;
    }

    /*
     * Period in minutes
     */
    double Period() const
    {
        return period_;
    }

    /*
     * EPOCH
     */
    Julian Epoch() const
    {
        return epoch_;
    }

private:

    double mean_anomoly_;
    double ascending_node_;
    double argument_perigee_;
    double eccentricity_;
    double inclination_;
    double mean_motion_;
    double bstar_;
    double recovered_semi_major_axis_;
    double recovered_mean_motion_;
    double perigee_;
    double period_;
    Julian epoch_;
};

#endif

