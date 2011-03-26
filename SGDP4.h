#ifndef SGDP4_H_
#define SGDP4_H_

#include "Tle.h"

class SGDP4 {
public:
    SGDP4(void);
    virtual ~SGDP4(void);

    void SetTle(const Tle& tle);

    void FindPosition(double tsince);

    struct TleData {
        double bstar;
        double eo;
        double omega;
        double xincl;
        double xmo;
        double xno;
        double xnodeo;
        Julian epoch;
    };

private:
    void Initialize();

    bool first_run_;

    struct TleData tle_data_0_;

    double cosio_;
    double sinio_;
    double aodp_;
    double x3thm1_;

    double eta_;
    double coef1_;
    double c1_;
    double a3ovk2_;
    double x1mth2_;
    double c4_;
    double c5_;

    double xmdot_;
    double omgdot_;
    double xnodot_;
    double xnodcf_;
    double t2cof_;
    double xlcof_;
    double aycof_;
    double x7thm1_;
    double omgcof_;
    double xmcof_;
    double delmo_;
    double sinmo_;

    double d2_;
    double d3_;
    double d4_;
    double t3cof_;
    double t4cof_;
    double t5cof_;

    double gsto_;

    bool use_simple_model_;
    bool use_deep_space_;

    /*
     * XMO
     */
    double MeanAnomoly() const {
        return mean_anomoly_;
    }

    /*
     * XNODEO
     */
    double AscendingNode() const {
        return ascending_node_;
    }

    /*
     * OMEGAO
     */
    double ArgumentPerigee() const {
        return argument_perigee_;
    }

    /*
     * EO
     */
    double Eccentricity() const {
        return eccentricity_;
    }

    /*
     * XINCL
     */
    double Inclination() const {
        return inclination_;
    }

    /*
     * XNO
     */
    double MeanMotion() const {
        return mean_motion_;
    }

    /*
     * BSTAR
     */
    double BStar() const {
        return bstar_;
    }

    /*
     * AODP
     */
    double RecoveredSemiMajorAxis() const {
        return recovered_semi_major_axis;
    }

    /*
     * XNODP
     */
    double RecoveredMeanMotion() const {
        return recovered_mean_motion_;
    }

    /*
     * EPOCH
     */
    Julian Epoch() const {
        return epoch_;
    }

    double mean_anomoly_;
    double ascending_node_;
    double argument_perigee_;
    double eccentricity_;
    double inclination_;
    double mean_motion_;
    double bstar_;
    double recovered_semi_major_axis;
    double recovered_mean_motion_;
    Julian epoch_;

};

#endif

