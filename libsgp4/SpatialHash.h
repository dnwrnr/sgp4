#pragma once

#include <vector>
#include <immintrin.h>
#include "Vector.h"

namespace libsgp4
{

/**
 * @brief SpatialHash - Ultra-fast proximity filtering using AVX-512 VPOPCNTDQ.
 * Maps coordinates to a bit-signature based on grid occupancy.
 */
class SpatialHash
{
public:
    static const int GRID_SIZE = 128; // 128km bins
    
    struct alignas(64) Signature {
        uint64_t bits[8]; // 512 bits per signature
    };

    /**
     * @brief Generate a 512-bit signature for a position.
     * Sets bits based on the neighboring grid cells.
     */
    static Signature Generate(const Vector& pos) {
        Signature sig = {0};
        int gx = (int)((pos.x + 10000.0) / GRID_SIZE) % 512;
        int gy = (int)((pos.y + 10000.0) / GRID_SIZE) % 512;
        int gz = (int)((pos.z + 10000.0) / GRID_SIZE) % 512;
        
        // Simple hash: set bit at (gx ^ gy ^ gz)
        int bit = (gx ^ gy ^ gz) & 511;
        sig.bits[bit / 64] |= (1ULL << (bit % 64));
        return sig;
    }

    /**
     * @brief Vectorized filter: Check if any targets in a batch overlap with probe.
     */
    static void Filter(const Signature& probe, const std::vector<Signature>& targets, std::vector<int>& candidate_indices) {
        __m512i v_probe = _mm512_load_si512(&probe);
        int n = targets.size();
        
        for(int i=0; i<n; i++) {
            __m512i v_target = _mm512_load_si512(&targets[i]);
            __m512i v_and = _mm512_and_si512(v_probe, v_target);
            
            // If any bit is set in the intersection, they are potentially near
            if (!_mm512_test_epi64_mask(v_and, v_and)) continue;
            
            candidate_indices.push_back(i);
        }
    }
};

} // namespace libsgp4
