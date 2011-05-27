#ifndef SGDP4_H_
#define SGDP4_H_

#include "Tle.h"
#include "Eci.h"

class SGP4 {
public:
    SGP4(void);
    SGP4(const Tle& tle);
    virtual ~SGP4(void);

    void SetTle(const Tle& tle);
    void FindPosition(Eci* eci, double tsince) const;
    void FindPosition(Eci* eci, const Julian& date) const;

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

    struct CommonConstants {
        double cosio;
        double sinio;
        double eta;
        double t2cof;
        double a3ovk2;
        double x1mth2;
        double x3thm1;
        double x7thm1;
        double aycof;
        double xlcof;
        double xnodcf;
        double c1;
        double c4;
        double omgdot; // secular rate of omega (radians/sec)
        double xnodot; // secular rate of xnode (radians/sec)
        double xmdot; // secular rate of xmo (radians/sec)
    };

    struct NearSpaceConstants {
        double c5;
        double omgcof;
        double xmcof;
        double delmo;
        double sinmo;
        double d2;
        double d3;
        double d4;
        double t3cof;
        double t4cof;
        double t5cof;
    };

    struct DeepSpaceConstants {
        double gsto;
        double zmol;
        double zmos;
        /*
         * whether the deep space orbit is
         * geopotential resonance for 12 hour orbits
         */
        bool resonance_flag;
        /*
         * whether the deep space orbit is
         * 24h synchronous resonance
         */
        bool synchronous_flag;
        /*
         * lunar / solar constants for epoch
         * applied during DeepSpaceSecular()
         */
        double sse;
        double ssi;
        double ssl;
        double ssg;
        double ssh;
        /*
         * lunar / solar constants
         * used during DeepSpaceCalculateLunarSolarTerms()
         */
        double se2;
        double si2;
        double sl2;
        double sgh2;
        double sh2;
        double se3;
        double si3;
        double sl3;
        double sgh3;
        double sh3;
        double sl4;
        double sgh4;
        double ee2;
        double e3;
        double xi2;
        double xi3;
        double xl2;
        double xl3;
        double xl4;
        double xgh2;
        double xgh3;
        double xgh4;
        double xh2;
        double xh3;
        /*
         * used during DeepSpaceCalcDotTerms()
         */
        double d2201;
        double d2211;
        double d3210;
        double d3222;
        double d4410;
        double d4422;
        double d5220;
        double d5232;
        double d5421;
        double d5433;
        double del1;
        double del2;
        double del3;
    };

    struct IntegratorConstants {
        /*
         * integrator constants
         */
        double xfact;
        double xlamo;

        /*
         * integrator values for epoch
         */
        double xndot_0;
        double xnddt_0;
        double xldot_0;
    };

    struct IntegratorParams {
        /*
         * integrator values
         */
        double xli;
        double xni;
        double atime;
        /*
         * itegrator values for current d_atime_
         */
        double xndot_t;
        double xnddt_t;
        double xldot_t;
    };

private:
    void Initialize();
    void DeepSpaceInitialize(const double& eosq, const double& sinio, const double& cosio, const double& betao,
            const double& theta2, const double& betao2,
            const double& xmdot, const double& omgdot, const double& xnodot);
    void DeepSpaceCalculateLunarSolarTerms(const double t, double* pe, double* pinc,
            double* pl, double* pgh, double* ph) const;
    void DeepSpacePeriodics(const double& t, double* em, double* xinc,
            double* omgasm, double* xnodes, double* xll) const;
    void DeepSpaceSecular(const double& t, double* xll, double* omgasm,
            double* xnodes, double* em, double* xinc, double* xn) const;
    void FindPositionSDP4(Eci* eci, double tsince) const;
    void FindPositionSGP4(Eci* eci, double tsince) const;
    void CalculateFinalPositionVelocity(Eci* eci, const double& tsince, const double& e,
            const double& a, const double& omega, const double& xl, const double& xnode,
            const double& xincl, const double& xlcof, const double& aycof,
            const double& x3thm1, const double& x1mth2, const double& x7thm1,
            const double& cosio, const double& sinio) const;
    void DeepSpaceCalcDotTerms(double* xndot, double* xnddt, double* xldot) const;
    void DeepSpaceIntegrator(const double delt, const double step2,
            const double xndot, const double xnddt, const double xldot)const;
    void Reset();

    bool first_run_;
    bool use_simple_model_;
    bool use_deep_space_;

    struct CommonConstants common_consts_;
    struct NearSpaceConstants nearspace_consts_;
    struct DeepSpaceConstants deepspace_consts_;
    struct IntegratorConstants integrator_consts_;
    mutable struct IntegratorParams integrator_params_;

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
};

#endif

