#pragma once

#include "Tle.h"
#include <vector>
#include "Eci.h"
#include "DateTime.h"
#include <memory>

namespace libsgp4
{

/**
 * @brief SGP4Batch - High-performance vectorized SGP4 propagator.
 */
class SGP4Batch
{
public:
    enum class MathMode {
        Standard,   // libmvec (High precision)
        Minimax     // Fast Polynomial Approx (~1e-12 precision)
    };

    SGP4Batch(const std::vector<Tle>& tles);
    void Propagate(double tsince, std::vector<Eci>& results, MathMode mode = MathMode::Standard) const;

    int total_satellites() const { return total_satellites_; }

    // JIT Accessors
    const double* get_cosio() const { return cosio.get(); }
    const double* get_sinio() const { return sinio.get(); }
    const double* get_eta() const { return eta.get(); }
    const double* get_t2cof() const { return t2cof.get(); }
    const double* get_x1mth2() const { return x1mth2.get(); }
    const double* get_x3thm1() const { return x3thm1.get(); }
    const double* get_x7thm1() const { return x7thm1.get(); }
    const double* get_aycof() const { return aycof.get(); }
    const double* get_xlcof() const { return xlcof.get(); }
    const double* get_xnodcf() const { return xnodcf.get(); }
    const double* get_c1() const { return c1.get(); }
    const double* get_c4() const { return c4.get(); }
    const double* get_omgdot() const { return omgdot.get(); }
    const double* get_xnodot() const { return xnodot.get(); }
    const double* get_xmdot() const { return xmdot.get(); }
    const double* get_c5() const { return c5.get(); }
    const double* get_omgcof() const { return omgcof.get(); }
    const double* get_xmcof() const { return xmcof.get(); }
    const double* get_delmo() const { return delmo.get(); }
    const double* get_sinmo() const { return sinmo.get(); }
    const double* get_d2() const { return d2.get(); }
    const double* get_d3() const { return d3.get(); }
    const double* get_d4() const { return d4.get(); }
    const double* get_t3cof() const { return t3cof.get(); }
    const double* get_t4cof() const { return t4cof.get(); }
    const double* get_t5cof() const { return t5cof.get(); }
    const double* get_xmo() const { return xmo.get(); }
    const double* get_nodeo() const { return nodeo.get(); }
    const double* get_omegao() const { return omegao.get(); }
    const double* get_ecco() const { return ecco.get(); }
    const double* get_inclo() const { return inclo.get(); }
    const double* get_bstar() const { return bstar.get(); }
    const double* get_aodp() const { return aodp.get(); }
    const double* get_no_kozai() const { return no_kozai.get(); }
    const double* get_active() const { return active.get(); }
    const std::vector<DateTime>& get_epochs() const { return epochs; }

private:
    int total_satellites_;
    int padded_satellites_;

    template<typename T>
    struct AlignedDeleter {
        void operator()(T* p) const { free(p); }
    };

    template<typename T>
    using AlignedPtr = std::unique_ptr<T[], AlignedDeleter<T>>;

    template<typename T>
    static AlignedPtr<T> make_aligned(size_t size) {
        void* ptr = nullptr;
        if (posix_memalign(&ptr, 64, size * sizeof(T)) != 0) throw std::bad_alloc();
        return AlignedPtr<T>(static_cast<T*>(ptr));
    }

    // --- Near Space SoA ---
    AlignedPtr<double> cosio, sinio, eta, t2cof, x1mth2, x3thm1, x7thm1, aycof, xlcof, xnodcf, c1, c4, omgdot, xnodot, xmdot;
    AlignedPtr<double> c5, omgcof, xmcof, delmo, sinmo, d2, d3, d4, t3cof, t4cof, t5cof;
    AlignedPtr<double> xmo, nodeo, omegao, ecco, inclo, bstar, aodp, no_kozai;
    AlignedPtr<double> use_simple_model, use_deep_space, active;

    // --- Deep Space SoA ---
    AlignedPtr<double> ds_gsto, ds_zmol, ds_zmos;
    AlignedPtr<double> ds_sse, ds_ssi, ds_ssl, ds_ssh, ds_ssg;
    AlignedPtr<double> ds_se2, ds_si2, ds_sl2, ds_sgh2, ds_sh2;
    AlignedPtr<double> ds_se3, ds_si3, ds_sl3, ds_sgh3, ds_sh3;
    AlignedPtr<double> ds_sl4, ds_sgh4, ds_ee2, ds_e3, ds_xi2, ds_xi3, ds_xl2, ds_xl3, ds_xl4, ds_xgh2, ds_xgh3, ds_xgh4, ds_xh2, ds_xh3;
    AlignedPtr<double> ds_d2201, ds_d2211, ds_d3210, ds_d3222, ds_d4410, ds_d4422, ds_d5220, ds_d5232, ds_d5421, ds_d5433;
    AlignedPtr<double> ds_del1, ds_del2, ds_del3, ds_xfact, ds_xlamo, ds_shape;

    // Integrator state
    mutable AlignedPtr<double> int_xli, int_xni, int_atime;

    std::vector<DateTime> epochs;

    // Memoization
    mutable AlignedPtr<double> last_tsince, memo_e, memo_a, memo_omega, memo_xl, memo_xnode;

    friend class JitPropagator;
};

} // namespace libsgp4
