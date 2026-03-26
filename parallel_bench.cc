#include <SGP4Batch.h>
#include <Tle.h>
#include <Vector.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <set>
#include <map>
#include <omp.h>
#include <chrono>

using namespace libsgp4;

struct Threat {
    std::string id;
    Tle tle;
    double peri;
    double apo;
};

inline double get_dist(const Vector& v1, const Vector& v2) {
    double dx = v1.x - v2.x, dy = v1.y - v2.y, dz = v1.z - v2.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

int main() {
    std::set<std::string> inactive_ids;
    std::ifstream ids_file("/home/troyrock/.openclaw/workspace/inactive_ids.txt");
    std::string id_line;
    while (std::getline(ids_file, id_line)) {
        id_line.erase(id_line.find_last_not_of(" \n\r\t")+1);
        if (!id_line.empty()) inactive_ids.insert(id_line);
    }
    std::cout << "Loaded " << inactive_ids.size() << " inactive IDs." << std::endl;

    std::map<std::string, std::pair<std::string, std::string>> latest_tles;
    std::ifstream tle_file("tle_sample.txt");
    std::string line, l1;
    while (std::getline(tle_file, line)) {
        if (line.empty()) continue;
        if (line[0] == '1') {
            l1 = line;
        } else if (line[0] == '2') {
            if (l1.empty()) continue;
            if (line.size() < 7) continue;
            std::string id = line.substr(2, 5);
            if (inactive_ids.count(id)) {
                latest_tles[id] = {l1, line};
            }
            l1.clear();
        }
    }
    
    std::vector<Tle> tles;
    std::vector<std::pair<double, double>> shell_bounds;
    for (auto const& [id, lines] : latest_tles) {
        try {
            std::string s1 = lines.first.substr(0, 69);
            std::string s2 = lines.second.substr(0, 69);
            Tle tle(id, s1, s2);
            double n = tle.MeanMotion();
            double a = std::pow(398600.44 / (n * n * (2.0*M_PI/86400.0) * (2.0*M_PI/86400.0)), 1.0/3.0);
            double e = tle.Eccentricity();
            tles.push_back(tle);
            shell_bounds.push_back({a*(1.0-e)-6378.137, a*(1.0+e)-6378.137});
        } catch (...) {}
    }
    std::cout << "Initialized " << tles.size() << " inactive threats." << std::endl;

    SGP4Batch batch_model(tles);
    std::cout << "Starting parallelized propagation benchmark..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    
    // We'll simulate 60 minutes of propagation for the entire batch
    // RUN 1: FRESH CALC
    std::vector<Eci> results;
    for (int t = 0; t < 60; ++t) {
        batch_model.Propagate(static_cast<double>(t), results);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    
    std::cout << "Fresh Run Throughput: " << (static_cast<double>(tles.size()) * 60.0 / diff.count() / 1e6) << " million sats/sec" << std::endl;

    // RUN 2: MEMOIZED
    start = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < 60; ++t) {
        batch_model.Propagate(static_cast<double>(t), results);
    }
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Memoized Run Throughput: " << (static_cast<double>(tles.size()) * 60.0 / diff.count() / 1e6) << " million sats/sec" << std::endl;

    return 0;
}
