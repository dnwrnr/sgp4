#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <algorithm>
#include <map>
#include "DateTime.h"
#include "Vector.h"
#include "Eci.h"

namespace libsgp4 {

/**
 * @brief TemporalPruner - Manages "sleep" states for satellite pairs based on 
 * minimum time-to-collision (TTC).
 */
class TemporalPruner {
public:
    struct PairState {
        uint32_t id1;
        uint32_t id2;
        int sleep_steps;
    };

    /**
     * @brief Estimates steps to skip based on distance and max closing velocity.
     * @param dist Current distance (km)
     * @param threshold_km Screening threshold (km)
     * @param step_size_sec Current simulation step (sec)
     * @param max_rel_vel_kms Conservative max relative velocity (typically 15.0 km/s in LEO)
     * @return Number of steps to safely sleep
     */
    static inline int EstimateSleepSteps(double dist, double threshold_km, double step_size_sec, double max_rel_vel_kms = 16.0) {
        if (dist <= threshold_km) return 0;
        
        // Time to reach threshold at max velocity: t = (dist - threshold) / v_max
        double safe_time_sec = (dist - threshold_km) / max_rel_vel_kms;
        int steps = static_cast<int>(safe_time_sec / step_size_sec);
        
        // Return steps - 1 for safety margin, capped to avoid long-term stale states
        return std::min(steps, 200); 
    }

    /**
     * @brief Fast pair-indexing for sleep map
     */
    static inline uint64_t GetPairKey(int id1, int id2) {
        if (id1 > id2) std::swap(id1, id2);
        return (static_cast<uint64_t>(id1) << 32) | static_cast<uint64_t>(id2);
    }
};

} // namespace libsgp4
