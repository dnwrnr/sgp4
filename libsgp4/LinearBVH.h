#pragma once

#include <vector>
#include <algorithm>
#include <immintrin.h>
#include <cstdint>
#include "Vector.h"

namespace libsgp4
{

/**
 * @brief LinearBVH - High-performance vectorized Bounding Volume Hierarchy.
 * Built on Morton Codes (Z-Order Curve) for O(N log N) build and O(log N) query.
 */
class LinearBVH
{
public:
    struct alignas(64) Node {
        float min_x, min_y, min_z, pad1;
        float max_x, max_y, max_z, pad2;
        int left_idx;   // If < 0, it's a leaf. Index to objects is ~left_idx
        int right_idx;  // If leaf, this is count or unused.
        int parent_idx;
        int pad3;
    };

    struct Object {
        int id;
        float x, y, z;
        uint32_t morton;
    };

    /**
     * @brief Build the LBVH from a set of points.
     */
    void Build(const std::vector<Object>& objects);

    /**
     * @brief Vectorized query: Find all objects within radius of center.
     */
    void Query(const Vector& center, double radius, std::vector<int>& results) const;

    static uint32_t GenerateMorton(float x, float y, float z);

private:
    std::vector<Node> nodes_;
    std::vector<Object> sorted_objects_;

    static uint32_t ExpandBits(uint32_t v) {
        v = (v * 0x00010001u) & 0xFF0000FFu;
        v = (v * 0x00000101u) & 0x0F00F00Fu;
        v = (v * 0x00000011u) & 0xC30C30C3u;
        v = (v * 0x00000005u) & 0x49249249u;
        return v;
    }

    int LongestCommonPrefix(int i, int j) const {
        if (j < 0 || j >= (int)sorted_objects_.size()) return -1;
        uint32_t a = sorted_objects_[i].morton;
        uint32_t b = sorted_objects_[j].morton;
        if (a == b) return 32 + __builtin_clz((uint32_t)(i ^ j));
        return __builtin_clz(a ^ b);
    }
};

} // namespace libsgp4
