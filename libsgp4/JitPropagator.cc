#include "JitPropagator.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <dlfcn.h>
#include <unistd.h>
#include <sys/stat.h>

namespace libsgp4
{

JitPropagator::JitPropagator(const SGP4Batch& batch) {
    std::string source = GenerateSource(batch);
    CompileAndLoad(source);
}

JitPropagator::~JitPropagator() {
    if (handle_) dlclose(handle_);
    if (!lib_path_.empty()) unlink(lib_path_.c_str());
}

std::string JitPropagator::GenerateSource(const SGP4Batch& batch) {
    std::stringstream ss;
    ss << "#include <immintrin.h>\n"
       << "#include <cmath>\n"
       << "#include <vector>\n"
       << "#include \"SGP4Batch.h\"\n"
       << "#include \"Globals.h\"\n"
       << "#include \"Eci.h\"\n"
       << "#include \"Vector.h\"\n"
       << "#include \"DateTime.h\"\n\n"
       << "using namespace libsgp4;\n\n"
       << "extern \"C\" {\n"
       << "__m512d _ZGVeN8v_sin(__m512d x);\n"
       << "__m512d _ZGVeN8v_cos(__m512d x);\n"
       << "__m512d _ZGVeN8vv_atan2(__m512d y, __m512d x);\n"
       << "}\n\n"
       << "struct Vec8 {\n"
       << "    __m512d v;\n"
       << "    Vec8() : v(_mm512_setzero_pd()) {}\n"
       << "    Vec8(__m512d val) : v(val) {}\n"
       << "    Vec8(double s) : v(_mm512_set1_pd(s)) {}\n"
       << "    static Vec8 load(const double* p) { return _mm512_load_pd(p); }\n"
       << "    void store(double* p) const { _mm512_store_pd(p, v); }\n"
       << "    Vec8 operator+(Vec8 b) const { return _mm512_add_pd(v, b.v); }\n"
       << "    Vec8 operator-(Vec8 b) const { return _mm512_sub_pd(v, b.v); }\n"
       << "    Vec8 operator*(Vec8 b) const { return _mm512_mul_pd(v, b.v); }\n"
       << "    Vec8 operator/(Vec8 b) const { return _mm512_div_pd(v, b.v); }\n"
       << "    Vec8 fmadd(Vec8 b, Vec8 c) const { return _mm512_fmadd_pd(v, b.v, c.v); }\n"
       << "};\n"
       << "static inline void v_sincos_vector(Vec8 x, Vec8& s, Vec8& c) { s.v = _ZGVeN8v_sin(x.v); c.v = _ZGVeN8v_cos(x.v); }\n"
       << "static inline Vec8 v_atan2_vector(Vec8 y, Vec8 x) { return Vec8(_ZGVeN8vv_atan2(y.v, x.v)); }\n"
       << "static inline Vec8 v_wrap_twopi(Vec8 x) {\n"
       << "    alignas(64) double a[8], r[8]; x.store(a);\n"
       << "    for(int i=0; i<8; i++) { r[i] = fmod(a[i], kTWOPI); if(r[i] < 0) r[i] += kTWOPI; }\n"
       << "    return Vec8::load(r);\n"
       << "}\n\n"
       << "extern \"C\" void sgp4_jit_kernel(double tsince, const SGP4Batch& batch, std::vector<Eci>& results) {\n"
       << "    const int n_pad = " << ((batch.total_satellites() + 7) & ~7) << ";\n"
       << "    const Vec8 v_tsince(tsince); const Vec8 v_tsq = v_tsince * v_tsince;\n"
       << "    const Vec8 v_tcube = v_tsq * v_tsince; const Vec8 v_tfour = v_tsince * v_tcube;\n"
       << "    const Vec8 v_one(1.0); const Vec8 v_zero(0.0);\n\n"
       << "    const double* __restrict__ p_xmdot = batch.get_xmdot(); const double* __restrict__ p_xmo = batch.get_xmo();\n"
       << "    const double* __restrict__ p_omgdot = batch.get_omgdot(); const double* __restrict__ p_omegao = batch.get_omegao();\n"
       << "    const double* __restrict__ p_xnodot = batch.get_xnodot(); const double* __restrict__ p_nodeo = batch.get_nodeo();\n"
       << "    const double* __restrict__ p_xnodcf = batch.get_xnodcf(); const double* __restrict__ p_c1 = batch.get_c1();\n"
       << "    const double* __restrict__ p_bstar = batch.get_bstar(); const double* __restrict__ p_c4 = batch.get_c4();\n"
       << "    const double* __restrict__ p_t2cof = batch.get_t2cof(); const double* __restrict__ p_eta = batch.get_eta();\n"
       << "    const double* __restrict__ p_xmcof = batch.get_xmcof(); const double* __restrict__ p_delmo = batch.get_delmo();\n"
       << "    const double* __restrict__ p_omgcof = batch.get_omgcof(); const double* __restrict__ p_d2 = batch.get_d2();\n"
       << "    const double* __restrict__ p_d3 = batch.get_d3(); const double* __restrict__ p_d4 = batch.get_d4();\n"
       << "    const double* __restrict__ p_c5 = batch.get_c5(); const double* __restrict__ p_sinmo = batch.get_sinmo();\n"
       << "    const double* __restrict__ p_t3cof = batch.get_t3cof(); const double* __restrict__ p_t4cof = batch.get_t4cof();\n"
       << "    const double* __restrict__ p_t5cof = batch.get_t5cof(); const double* __restrict__ p_aodp = batch.get_aodp();\n"
       << "    const double* __restrict__ p_ecco = batch.get_ecco(); const double* __restrict__ p_no_kozai = batch.get_no_kozai();\n"
       << "    const double* __restrict__ p_xlcof = batch.get_xlcof(); const double* __restrict__ p_aycof = batch.get_aycof();\n"
       << "    const double* __restrict__ p_x3thm1 = batch.get_x3thm1(); const double* __restrict__ p_cosio = batch.get_cosio();\n"
       << "    const double* __restrict__ p_inclo = batch.get_inclo(); const double* __restrict__ p_sinio = batch.get_sinio();\n"
       << "    const double* __restrict__ p_x1mth2 = batch.get_x1mth2(); const double* __restrict__ p_x7thm1 = batch.get_x7thm1();\n"
       << "    const double* __restrict__ p_active = batch.get_active();\n\n"
       << "    #pragma omp parallel for schedule(static)\n"
       << "    for (int i = 0; i < n_pad; i += 8) {\n"
       << "        Vec8 xmdf = Vec8::load(&p_xmdot[i]).fmadd(v_tsince, Vec8::load(&p_xmo[i]));\n"
       << "        Vec8 omgadf = Vec8::load(&p_omgdot[i]).fmadd(v_tsince, Vec8::load(&p_omegao[i]));\n"
       << "        Vec8 xnoddf = Vec8::load(&p_xnodot[i]).fmadd(v_tsince, Vec8::load(&p_nodeo[i]));\n"
       << "        Vec8 xnode = xnoddf + Vec8::load(&p_xnodcf[i]) * v_tsq;\n"
       << "        Vec8 tempa = v_one - Vec8::load(&p_c1[i]) * v_tsince;\n"
       << "        Vec8 tempe = Vec8::load(&p_bstar[i]) * Vec8::load(&p_c4[i]) * v_tsince;\n"
       << "        Vec8 templ = Vec8::load(&p_t2cof[i]) * v_tsq;\n\n"
       << "        Vec8 sw_xmdf, cw_xmdf; v_sincos_vector(xmdf, sw_xmdf, cw_xmdf);\n"
       << "        Vec8 etacos = v_one + Vec8::load(&p_eta[i]) * cw_xmdf;\n"
       << "        Vec8 delm = Vec8::load(&p_xmcof[i]) * (etacos * etacos * etacos - Vec8::load(&p_delmo[i]));\n"
       << "        Vec8 delomg = Vec8::load(&p_omgcof[i]) * v_tsince;\n"
       << "        Vec8 xmp = xmdf + delomg + delm; Vec8 omega = omgadf - delomg - delm;\n"
       << "        tempa = tempa - Vec8::load(&p_d2[i])*v_tsq - Vec8::load(&p_d3[i])*v_tcube - Vec8::load(&p_d4[i])*v_tfour;\n"
       << "        Vec8 sw_xmp, cw_xmp; v_sincos_vector(xmp, sw_xmp, cw_xmp);\n"
       << "        tempe = tempe + Vec8::load(&p_bstar[i]) * Vec8::load(&p_c5[i]) * (sw_xmp - Vec8::load(&p_sinmo[i]));\n"
       << "        templ = templ + Vec8::load(&p_t3cof[i])*v_tcube + v_tfour * Vec8::load(&p_t5cof[i]).fmadd(v_tsince, Vec8::load(&p_t4cof[i]));\n\n"
       << "        Vec8 a = Vec8::load(&p_aodp[i]) * tempa * tempa;\n"
       << "        Vec8 e = Vec8::load(&p_ecco[i]) - tempe;\n"
       << "        e.v = _mm512_max_pd(e.v, _mm512_set1_pd(1e-6)); e.v = _mm512_min_pd(e.v, _mm512_set1_pd(1.0 - 1e-6));\n"
       << "        Vec8 xl = xmp + omega + xnode + Vec8::load(&p_no_kozai[i]) * templ;\n\n"
       << "        Vec8 beta2 = v_one - e * e; Vec8 sw_om, cw_om; v_sincos_vector(omega, sw_om, cw_om);\n"
       << "        Vec8 axn = e * cw_om; Vec8 temp11 = v_one / (a * beta2);\n"
       << "        Vec8 xlt = xl + temp11 * Vec8::load(&p_xlcof[i]) * axn;\n"
       << "        Vec8 ayn = e * sw_om + temp11 * Vec8::load(&p_aycof[i]);\n"
       << "        Vec8 capu = v_wrap_twopi(xlt - xnode); Vec8 epw = capu; Vec8 sinepw, cosepw, ecose, esine;\n"
       << "        for(int j=0; j<3; j++) {\n"
       << "            v_sincos_vector(epw, sinepw, cosepw);\n"
       << "            ecose = axn * cosepw + ayn * sinepw; esine = axn * sinepw - ayn * cosepw;\n"
       << "            epw = epw + (capu - epw + esine) / (v_one - ecose + Vec8(0.5) * esine * ((capu - epw + esine) / (v_one - ecose)));\n"
       << "        }\n"
       << "        v_sincos_vector(epw, sinepw, cosepw); ecose = axn * cosepw + ayn * sinepw; esine = axn * sinepw - ayn * cosepw;\n"
       << "        Vec8 r = a * (v_one - ecose); Vec8 temp31 = v_one / r;\n"
       << "        Vec8 rdot = Vec8(kXKE) * _mm512_sqrt_pd(a.v) * esine * temp31;\n"
       << "        Vec8 rfdot = Vec8(kXKE) * _mm512_sqrt_pd((a * (v_one - axn*axn - ayn*ayn)).v) * temp31;\n"
       << "        Vec8 cosu_un = a * temp31 * (cosepw - axn + ayn * esine * (v_one / (v_one + _mm512_sqrt_pd((v_one - axn*axn - ayn*ayn).v))));\n"
       << "        Vec8 sinu_un = a * temp31 * (sinepw - ayn - axn * esine * (v_one / (v_one + _mm512_sqrt_pd((v_one - axn*axn - ayn*ayn).v))));\n"
       << "        Vec8 u = v_atan2_vector(sinu_un, cosu_un);\n"
       << "        Vec8 sin2u = Vec8(2.0) * sinu_un * cosu_un; Vec8 cos2u = Vec8(2.0) * cosu_un * cosu_un - v_one;\n"
       << "        Vec8 t41 = v_one / (a * (v_one - axn*axn - ayn*ayn)); Vec8 t42 = Vec8(kCK2) * t41; Vec8 t43 = t42 * t41;\n"
       << "        Vec8 rk = r * (v_one - Vec8(1.5) * t43 * _mm512_sqrt_pd((v_one - axn*axn - ayn*ayn).v) * Vec8::load(&p_x3thm1[i])) + Vec8(0.5) * t42 * Vec8::load(&p_x1mth2[i]) * cos2u;\n"
       << "        Vec8 uk = u - Vec8(0.25) * t43 * Vec8::load(&p_x7thm1[i]) * sin2u;\n"
       << "        Vec8 xnodek = xnode + Vec8(1.5) * t43 * Vec8::load(&p_cosio[i]) * sin2u;\n"
       << "        Vec8 xinck = Vec8::load(&p_inclo[i]) + Vec8(1.5) * t43 * Vec8::load(&p_cosio[i]) * Vec8::load(&p_sinio[i]) * cos2u;\n"
       << "        Vec8 rdotk = rdot - Vec8(kXKE) / (a * _mm512_sqrt_pd(a.v)) * t42 * Vec8::load(&p_x1mth2[i]) * sin2u;\n"
       << "        Vec8 rfdotk = rfdot + Vec8(kXKE) / (a * _mm512_sqrt_pd(a.v)) * t42 * (Vec8::load(&p_x1mth2[i]) * cos2u + Vec8(1.5) * Vec8::load(&p_x3thm1[i]));\n"
       << "        Vec8 snuk, csuk, snik, csik, snnk, csnk;\n"
       << "        v_sincos_vector(uk, snuk, csuk); v_sincos_vector(xinck, snik, csik); v_sincos_vector(xnodek, snnk, csnk);\n"
       << "        Vec8 xmx = (v_zero - snnk) * csik; Vec8 xmy = csnk * csik;\n"
       << "        Vec8 ux = xmx * snuk + csnk * csuk; Vec8 uy = xmy * snuk + snnk * csuk; Vec8 uz = snik * snuk;\n"
       << "        Vec8 vx = xmx * csuk - csnk * snuk; Vec8 vy = xmy * csuk - snnk * snuk; Vec8 vz = snik * csuk;\n"
       << "        alignas(64) double a_rk[8], a_ux[8], a_uy[8], a_uz[8], a_rdotk[8], a_rfdotk[8], a_vx[8], a_vy[8], a_vz[8];\n"
       << "        rk.store(a_rk); ux.store(a_ux); uy.store(a_uy); uz.store(a_uz); rdotk.store(a_rdotk); rfdotk.store(a_rfdotk); vx.store(a_vx); vy.store(a_vy); vz.store(a_vz);\n"
       << "        for (int s = 0; s < 8; s++) {\n"
       << "            int idx = i + s; if (idx >= " << batch.total_satellites() << ") break;\n"
       << "            if (p_active[idx] == 0.0) continue;\n"
       << "            results[idx] = Eci(batch.get_epochs()[idx].AddMinutes(tsince), Vector(a_rk[s] * a_ux[s] * kXKMPER, a_rk[s] * a_uy[s] * kXKMPER, a_rk[s] * a_uz[s] * kXKMPER),\n"
       << "                           Vector((a_rdotk[s] * a_ux[s] + a_rfdotk[s] * a_vx[s]) * kXKMPER / 60.0, (a_rdotk[s] * a_uy[s] + a_rfdotk[s] * a_vy[s]) * kXKMPER / 60.0, (a_rdotk[s] * a_uz[s] + a_rfdotk[s] * a_vz[s]) * kXKMPER / 60.0));\n"
       << "        }\n"
       << "    }\n"
       << "}\n";
    return ss.str();
}

void JitPropagator::CompileAndLoad(const std::string& source) {
    char tmp_cc[] = "/tmp/sgp4_jit_XXXXXX.cc";
    int fd = mkstemps(tmp_cc, 3);
    if (fd == -1) return;
    close(fd);

    std::ofstream ofs(tmp_cc);
    ofs << source;
    ofs.close();

    char tmp_so[] = "/tmp/sgp4_jit_XXXXXX.so";
    int fd_so = mkstemps(tmp_so, 3);
    if (fd_so == -1) return;
    close(fd_so);

    std::string cmd = "g++ -O3 -ffast-math -march=native -fopenmp -shared -fPIC -I. -I/home/troyrock/fastSGP4/libsgp4 " + std::string(tmp_cc) + " -o " + std::string(tmp_so) + " -lm";
    if (system(cmd.c_str()) != 0) {
        std::cerr << "JIT Compilation Failed" << std::endl;
        return;
    }

    handle_ = dlopen(tmp_so, RTLD_NOW);
    if (!handle_) {
        std::cerr << "dlopen failed: " << dlerror() << std::endl;
        return;
    }

    kernel_func_ = (KernelFunc)dlsym(handle_, "sgp4_jit_kernel");
    lib_path_ = tmp_so;
    unlink(tmp_cc);
}

void JitPropagator::Propagate(double tsince, const SGP4Batch& batch, std::vector<Eci>& results) {
    if (kernel_func_) {
        kernel_func_(tsince, batch, results);
    }
}

} // namespace libsgp4
