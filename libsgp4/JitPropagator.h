#pragma once

#include "SGP4Batch.h"
#include <string>
#include <vector>
#include <dlfcn.h>
#include <functional>

namespace libsgp4
{

/**
 * @brief JitPropagator - Runtime specialized kernel generator for SGP4Batch.
 */
class JitPropagator
{
public:
    using KernelFunc = void(*)(double tsince, const SGP4Batch& batch, std::vector<Eci>& results);

    JitPropagator(const SGP4Batch& batch);
    ~JitPropagator();

    /**
     * @brief Propagate using the JIT-compiled specialized kernel.
     */
    void Propagate(double tsince, const SGP4Batch& batch, std::vector<Eci>& results);

    /**
     * @brief Check if a specialized kernel is currently loaded and ready.
     */
    bool IsReady() const { return kernel_func_ != nullptr; }

private:
    std::string GenerateSource(const SGP4Batch& batch);
    void CompileAndLoad(const std::string& source);

    void* handle_ = nullptr;
    KernelFunc kernel_func_ = nullptr;
    std::string lib_path_;
};

} // namespace libsgp4
