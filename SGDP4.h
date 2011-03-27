#ifndef SGDP4_H_
#define SGDP4_H_

#include "Tle.h"

class SGDP4 {
public:
    SGDP4(void);
    virtual ~SGDP4(void);

    void SetTle(const Tle& tle);
    void FindPosition(double tsince);

private:
    void Initialize(const double& theta2, const double& betao2, const double& betao, const double& eosq);
    void DeepSpaceInitialize(const double& eosq, const double& sinio, const double& cosio, const double& betao,
            const double& theta2, const double& sing, const double& cosg, const double& betao2,
            const double& xmdot, const double& omgdot, const double& xnodot);
    void DeepPeriodics(const double& sinio, const double& cosio, const double& t, double& em, double& xinc,
            double& omgasm, double& xnodes, double& xll);

    bool first_run_;

    double cosio_;
    double sinio_;
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
        return recovered_semi_major_axis_;
    }

    /*
     * XNODP
     */
    double RecoveredMeanMotion() const {
        return recovered_mean_motion_;
    }

    /*
     * PERIGE
     */
    double Perigee() const {
        return perigee_;
    }

    double Period() const {
        return period_;
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
    double recovered_semi_major_axis_;
    double recovered_mean_motion_;
    double perigee_;
    double period_;
    Julian epoch_;














    double d_day_;

    double d_xnodce_;
    double d_stem_;
    double d_ctem_;
    double d_zcosil_;
    double d_zsinil_;
    double d_zsinhl_;
    double d_zcoshl_;
    double d_c_;
    double d_gam_;
    double d_zx_;
    double d_zy_;
    double d_zcosgl_;
    double d_zsingl_;
    double d_zmol_;
    double d_zmos_;

    double d_savtsn_;

    double d_sse_;
    double d_ssi_;
    double d_ssl_;
    double d_ssg_;
    double d_ssh_;
    double d_se2_;
    double d_si2_;
    double d_sl2_;
    double d_sgh2_;
    double d_sh2_;
    double d_se3_;
    double d_si3_;
    double d_sl3_;
    double d_sgh3_;
    double d_sh3_;
    double d_sl4_;
    double d_sgh4_;

    double d_ee2_;
    double d_e3_;
    double d_xi2_;
    double d_xi3_;
    double d_xl2_;
    double d_xl3_;
    double d_xl4_;
    double d_xgh2_;
    double d_xgh3_;
    double d_xgh4_;
    double d_xh2_;
    double d_xh3_;

    double d_d2201_;
    double d_d2211_;
    double d_d3210_;
    double d_d3222_;
    double d_d4410_;
    double d_d4422_;
    double d_d5220_;
    double d_d5232_;
    double d_d5421_;
    double d_d5433_;

    double d_del1_;
    double d_del2_;
    double d_del3_;
    double d_fasx2_;
    double d_fasx4_;
    double d_fasx6_;

    double d_xfact_;

    double d_xlamo_;
    double d_xli_;
    double d_xni_;
    double d_atime_;
    double d_stepp_;
    double d_stepn_;
    double d_step2_;

};

#endif

