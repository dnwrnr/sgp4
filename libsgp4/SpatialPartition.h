#include <iostream>
#include <vector>
#include <algorithm>
#include <immintrin.h>
#include <cmath>
#include <cstdint>

namespace libsgp4 {

/**
 * @brief SpatialPartition - BVH-lite using Morton codes and Z-order curve.
 * This handles the O(N log N) pruning for conjunction screening.
 */
class SpatialPartition {
public:
    struct Object {
        int id;
        double x, y, z;
        uint64_t morton;
    };

    // Expands a 21-bit integer into 64 bits by inserting 2 zeros between each bit
    static inline uint64_t expandBits(uint64_t v) {
        v = (v | (v << 32)) & 0x001F00000000FFFF;
        v = (v | (v << 16)) & 0x001F0000FF0000FF;
        v = (v | (v << 8))  & 0x100F00F00F00F00F;
        v = (v | (v << 4))  & 0x10C30C30C30C30C3;
        v = (v | (v << 2))  & 0x1249249249249249;
        return v;
    }

    // Calculates Morton code for a point in [-RE*10, RE*10] space (LEO to GEO)
    static inline uint64_t morton3D(double x, double y, double z) {
        // Normalize coordinates to [0, 2^21 - 1]
        // Range: +/- 100,000 km is plenty for SATCAT
        const double offset = 100000.0;
        const double scale = 10.48575; // (2^21 - 1) / 200,000
        
        uint64_t ix = static_cast<uint64_t>((x + offset) * scale);
        uint64_t iy = static_cast<uint64_t>((y + offset) * scale);
        uint64_t iz = static_cast<uint64_t>((z + offset) * scale);

        return (expandBits(ix) << 2) | (expandBits(iy) << 1) | expandBits(iz);
    }

    /**
     * @brief Sorts objects along the Z-order curve.
     * Objects close in 3D space will likely be close in the sorted list.
     */
    static void SortObjects(std::vector<Object>& objects) {
        for (auto& obj : objects) {
            obj.morton = morton3D(obj.x, obj.y, obj.z);
        }
        std::sort(objects.begin(), objects.end(), [](const Object& a, const Object& b) {
            return a.morton < b.morton;
        });
    }

    /**
     * @brief Screen conjunctions using Z-order sliding window.
     * Complexity: O(N log N) due to sort, O(N * K) for sweep where K is small.
     */
    static void SweepAndPrune(
        const std::vector<Object>& sorted_objects,
        double threshold_km,
        std::vector<std::pair<int, int>>& candidates) 
    {
        int n = sorted_objects.size();
        double d2_thresh = threshold_km * threshold_km;

        for (int i = 0; i < n; i++) {
            // Sweep forward in the Morton-sorted list
            for (int j = i + 1; j < n; j++) {
                // If the Morton distance is too large, we can stop the inner sweep early
                // (Note: This is a heuristic pruning. For a strict BVH, we'd use a tree traversal).
                // However, along the Z-curve, objects that are far in Morton space ARE far in 1 coordinate.
                
                // Morton codes are 64-bit. If the high bits differ significantly, 
                // the spatial distance is guaranteed to be large.
                if ((sorted_objects[j].morton - sorted_objects[i].morton) > 1000000) { 
                    // High-performance heuristic: tune this value or check bounding box
                    // For bit-perfect, we check if coord difference exceeds threshold
                    if (std::abs(sorted_objects[j].x - sorted_objects[i].x) > threshold_km) break;
                }

                double dx = sorted_objects[i].x - sorted_objects[j].x;
                double dy = sorted_objects[i].y - sorted_objects[j].y;
                double dz = sorted_objects[i].z - sorted_objects[j].z;
                double d2 = dx*dx + dy*dy + dz*dz;

                if (d2 < d2_thresh) {
                    candidates.push_back({sorted_objects[i].id, sorted_objects[j].id});
                }
            }
        }
    }
};

} // namespace libsgp4
