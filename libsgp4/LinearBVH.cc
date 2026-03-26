#include "LinearBVH.h"
#include <stack>
#include <immintrin.h>
#include <functional>
#include <cstdint>

namespace libsgp4
{

uint32_t LinearBVH::GenerateMorton(float x, float y, float z) {
    // Normalize to [0, 1023]
    uint32_t ux = (uint32_t)std::max(0.0f, std::min(1023.0f, (x + 15000.0f) * (1023.0f / 30000.0f)));
    uint32_t uy = (uint32_t)std::max(0.0f, std::min(1023.0f, (y + 15000.0f) * (1023.0f / 30000.0f)));
    uint32_t uz = (uint32_t)std::max(0.0f, std::min(1023.0f, (z + 15000.0f) * (1023.0f / 30000.0f)));
    return (ExpandBits(ux) << 2) | (ExpandBits(uy) << 1) | ExpandBits(uz);
}

void LinearBVH::Build(const std::vector<Object>& objects) {
    if (objects.empty()) return;
    int n = static_cast<int>(objects.size());
    sorted_objects_ = objects;
    
    // 1. Sort by Morton Code
    std::sort(sorted_objects_.begin(), sorted_objects_.end(), [](const Object& a, const Object& b){
        return a.morton < b.morton;
    });

    // 2. Build Hierarchy
    nodes_.assign(2 * n - 1, Node());
    
    std::function<int(int, int)> find_split = [&](int first, int last) {
        uint32_t first_code = sorted_objects_[static_cast<size_t>(first)].morton;
        uint32_t last_code = sorted_objects_[static_cast<size_t>(last)].morton;
        if (first_code == last_code) return (first + last) >> 1;

        int common_prefix = __builtin_clz(first_code ^ last_code);
        int split = first;
        int step = last - first;
        do {
            step = (step + 1) >> 1;
            int new_split = split + step;
            if (new_split < last) {
                int split_prefix = __builtin_clz(first_code ^ sorted_objects_[static_cast<size_t>(new_split)].morton);
                if (split_prefix > common_prefix) split = new_split;
            }
        } while (step > 1);
        return split;
    };

    std::function<int(int, int)> build_recursive = [&](int first, int last) {
        if (first == last) {
            int idx = n - 1 + first;
            nodes_[static_cast<size_t>(idx)].left_idx = ~first; 
            nodes_[static_cast<size_t>(idx)].min_x = sorted_objects_[static_cast<size_t>(first)].x; 
            nodes_[static_cast<size_t>(idx)].max_x = sorted_objects_[static_cast<size_t>(first)].x;
            nodes_[static_cast<size_t>(idx)].min_y = sorted_objects_[static_cast<size_t>(first)].y; 
            nodes_[static_cast<size_t>(idx)].max_y = sorted_objects_[static_cast<size_t>(first)].y;
            nodes_[static_cast<size_t>(idx)].min_z = sorted_objects_[static_cast<size_t>(first)].z; 
            nodes_[static_cast<size_t>(idx)].max_z = sorted_objects_[static_cast<size_t>(first)].z;
            return idx;
        }

        int split = find_split(first, last);
        int left = build_recursive(first, split);
        int right = build_recursive(split + 1, last);

        int idx = first; 
        nodes_[static_cast<size_t>(idx)].left_idx = left;
        nodes_[static_cast<size_t>(idx)].right_idx = right;
        nodes_[static_cast<size_t>(idx)].min_x = std::min(nodes_[static_cast<size_t>(left)].min_x, nodes_[static_cast<size_t>(right)].min_x);
        nodes_[static_cast<size_t>(idx)].max_x = std::max(nodes_[static_cast<size_t>(left)].max_x, nodes_[static_cast<size_t>(right)].max_x);
        nodes_[static_cast<size_t>(idx)].min_y = std::min(nodes_[static_cast<size_t>(left)].min_y, nodes_[static_cast<size_t>(right)].min_y);
        nodes_[static_cast<size_t>(idx)].max_y = std::max(nodes_[static_cast<size_t>(left)].max_y, nodes_[static_cast<size_t>(right)].max_y);
        nodes_[static_cast<size_t>(idx)].min_z = std::min(nodes_[static_cast<size_t>(left)].min_z, nodes_[static_cast<size_t>(right)].min_z);
        nodes_[static_cast<size_t>(idx)].max_z = std::max(nodes_[static_cast<size_t>(left)].max_z, nodes_[static_cast<size_t>(right)].max_z);
        return idx;
    };

    build_recursive(0, n - 1);
}

void LinearBVH::Query(const Vector& center, double radius, std::vector<int>& results) const {
    if (nodes_.empty()) return;
    
    float cx = (float)center.x, cy = (float)center.y, cz = (float)center.z;
    float r2 = (float)(radius * radius);

    int stack[64];
    int top = 0;
    stack[top++] = 0;

    while (top > 0) {
        int idx = stack[--top];
        const Node& node = nodes_[static_cast<size_t>(idx)];

        float dx = std::max(0.0f, std::max(node.min_x - cx, cx - node.max_x));
        float dy = std::max(0.0f, std::max(node.min_y - cy, cy - node.max_y));
        float dz = std::max(0.0f, std::max(node.min_z - cz, cz - node.max_z));
        
        if (dx*dx + dy*dy + dz*dz > r2) continue;

        if (node.left_idx < 0) { 
            int obj_idx = ~node.left_idx;
            results.push_back(sorted_objects_[static_cast<size_t>(obj_idx)].id);
        } else {
            stack[top++] = node.right_idx;
            stack[top++] = node.left_idx;
        }
    }
}

} // namespace libsgp4
