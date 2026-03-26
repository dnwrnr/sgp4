#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <immintrin.h>
#include "DateTime.h"
#include "Vector.h"
#include "Eci.h"

namespace libsgp4 {

/**
 * @brief HermiteInterpolator - Fast vectorized state interpolation.
 * Given states at T0 and T1, interpolates positions at T0 + delta.
 */
class HermiteInterpolator {
public:
    /**
     * @brief Interpolate a batch of satellites using Cubic Hermite Splines.
     * @param results0 States at start time (t0)
     * @param results1 States at end time (t1)
     * @param dt_total Total interval (t1 - t0) in seconds
     * @param dt_target Offset from t0 in seconds
     * @param out_positions Vector to store interpolated positions
     */
    static void InterpolatePositions(
        const std::vector<Eci>& results0,
        const std::vector<Eci>& results1,
        double dt_total,
        double dt_target,
        std::vector<Vector>& out_positions) 
    {
        int n = results0.size();
        out_positions.resize(n);

        double u = dt_target / dt_total;
        double u2 = u * u;
        double u3 = u2 * u;

        // Hermite Basis Functions
        double h00 = 2*u3 - 3*u2 + 1;
        double h10 = u3 - 2*u2 + u;
        double h01 = -2*u3 + 3*u2;
        double h11 = u3 - u2;

        // Multipliers for velocities (chain rule: dP/dt * dt/du)
        double v_mult0 = h10 * dt_total;
        double v_mult1 = h11 * dt_total;

        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            const Vector& p0 = results0[i].Position();
            const Vector& v0 = results0[i].Velocity(); // km/s
            const Vector& p1 = results1[i].Position();
            const Vector& v1 = results1[i].Velocity();

            out_positions[i] = Vector(
                h00 * p0.x + v_mult0 * v0.x + h01 * p1.x + v_mult1 * v1.x,
                h00 * p0.y + v_mult0 * v0.y + h01 * p1.y + v_mult1 * v1.y,
                h00 * p0.z + v_mult0 * v0.z + h01 * p1.z + v_mult1 * v1.z
            );
        }
    }

    /**
     * @brief AVX-512 version of the interpolator.
     * Requires positions/velocities in SoA layout.
     */
    static void InterpolateBatchSIMD(
        int n_sat,
        double dt_total,
        double dt_target,
        const double* p0x, const double* p0y, const double* p0z,
        const double* v0x, const double* v0y, const double* v0z,
        const double* p1x, const double* p1y, const double* p1z,
        const double* v1x, const double* v1y, const double* v1z,
        double* outx, double* outy, double* outz)
    {
        double u = dt_target / dt_total;
        double u2 = u * u;
        double u3 = u2 * u;

        __m512d v_h00 = _mm512_set1_pd(2*u3 - 3*u2 + 1);
        __m512d v_h10 = _mm512_set1_pd(u3 - 2*u2 + u);
        __m512d v_h01 = _mm512_set1_pd(-2*u3 + 3*u2);
        __m512d v_h11 = _mm512_set1_pd(u3 - u2);
        
        __m512d v_mult0 = _mm512_mul_pd(v_h10, _mm512_set1_pd(dt_total));
        __m512d v_mult1 = _mm512_mul_pd(v_h11, _mm512_set1_pd(dt_total));

        #pragma omp parallel for
        for (int i = 0; i < n_sat; i += 8) {
            // X
            __m512d vx = _mm512_mul_pd(v_h00, _mm512_loadu_pd(&p0x[i]));
            vx = _mm512_fmadd_pd(v_mult0, _mm512_loadu_pd(&v0x[i]), vx);
            vx = _mm512_fmadd_pd(v_h01, _mm512_loadu_pd(&p1x[i]), vx);
            vx = _mm512_fmadd_pd(v_mult1, _mm512_loadu_pd(&v1x[i]), vx);
            _mm512_storeu_pd(&outx[i], vx);

            // Y
            __m512d vy = _mm512_mul_pd(v_h00, _mm512_loadu_pd(&p0y[i]));
            vy = _mm512_fmadd_pd(v_mult0, _mm512_loadu_pd(&v0y[i]), vy);
            vy = _mm512_fmadd_pd(v_h01, _mm512_loadu_pd(&p1y[i]), vy);
            vy = _mm512_fmadd_pd(v_mult1, _mm512_loadu_pd(&v1y[i]), vy);
            _mm512_storeu_pd(&outy[i], vy);

            // Z
            __m512d vz = _mm512_mul_pd(v_h00, _mm512_loadu_pd(&p0z[i]));
            vz = _mm512_fmadd_pd(v_mult0, _mm512_loadu_pd(&v0z[i]), vz);
            vz = _mm512_fmadd_pd(v_h01, _mm512_loadu_pd(&p1z[i]), vz);
            vz = _mm512_fmadd_pd(v_mult1, _mm512_loadu_pd(&v1z[i]), vz);
            _mm512_storeu_pd(&outz[i], vz);
        }
    }
};

} // namespace libsgp4
