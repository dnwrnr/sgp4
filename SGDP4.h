#ifndef SGDP4_H_
#define SGDP4_H_

#include "Tle.h"
#include "Eci.h"

class SGDP4 {
public:
    SGDP4(void);
    virtual ~SGDP4(void);

    enum EnumConstants {
        CONSTANTS_WGS72_OLD,
        CONSTANTS_WGS72,
        CONSTANTS_WGS84
    };

    void SetConstants(EnumConstants constants);
    void SetTle(const Tle& tle);
    void FindPosition(Eci& eci, double tsince);

private:
    void Initialize(const double& theta2, const double& betao2, const double& betao, const double& eosq);
    void DeepSpaceInitialize(const double& eosq, const double& sinio, const double& cosio, const double& betao,
            const double& theta2, const double& sing, const double& cosg, const double& betao2,
            const double& xmdot, const double& omgdot, const double& xnodot);
    void DeepPeriodics(const double& t, double& em, double& xinc,
            double& omgasm, double& xnodes, double& xll);
    void DeepSecular(const double& t, double& xll, double& omgasm,
            double& xnodes, double& em, double& xinc, double& xn);
    void FindPositionSDP4(Eci& eci, double tsince);
    void FindPositionSGP4(Eci& eci, double tsince);
    void CalculateFinalPositionVelocity(Eci& eci, const double& tsince, const double& e,
            const double& a, const double& omega, const double& xl, const double& xnode,
            const double& xincl, const double& xlcof, const double& aycof,
            const double& x3thm1, const double& x1mth2, const double& x7thm1,
            const double& cosio, const double& sinio);
    void DeepSpaceCalcDotTerms(double& xndot, double& xnddt, double& xldot);
    void DeepSpaceCalcIntegrator(const double& delt, const double& step2, double& xndot, double& xnddt, double& xldot);

    bool first_run_;

    /*
     * variables are constants that wont be modified outside init
     */
    double i_cosio_;
    double i_sinio_;
    double i_eta_;
    double i_t2cof_;
    double i_a3ovk2_;
    double i_x1mth2_;
    double i_x3thm1_;
    double i_x7thm1_;
    double i_aycof_;
    double i_xlcof_;
    double i_xnodcf_;
    double i_c1_;
    double i_c4_;
    double i_omgdot_; // secular rate of omega (radians/sec)
    double i_xnodot_; // secular rate of xnode (radians/sec)
    double i_xmdot_; // secular rate of xmo (radians/sec)
    /*
     * sgp4 constant
     */
    double i_c5_;
    double i_omgcof_;
    double i_xmcof_;
    double i_delmo_;
    double i_sinmo_;
    double i_d2_;
    double i_d3_;
    double i_d4_;
    double i_t3cof_;
    double i_t4cof_;
    double i_t5cof_;
    /*
     * sdp4 constant
     */
    double i_gsto_;

    bool i_use_simple_model_;
    bool i_use_deep_space_;

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

    /*
     * these variables are set at the very start
     * and should not be changed after that
     */
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












    bool d_resonance_flag_;
    bool d_synchronous_flag_;

    double d_zmol_;
    double d_zmos_;

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

    struct Constants {
        double AE;
        double MU;
        double CK2;
        double CK4;
        double TWOTHRD;
        double XKE;
        double XKMPER;
        double S;
        double QOMS2T;
        double XJ2;
        double XJ3;
        double XJ4;
        double J3OJ2;
    };

    struct Constants constants_;

};

#endif

