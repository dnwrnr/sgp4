#include "libsgp4/SGP4.h"
#include "libsgp4/SGP4Batch.h"
#include "libsgp4/Tle.h"
#include "libsgp4/DateTime.h"
#include "libsgp4/Vector.h"
#include <chrono>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

using namespace libsgp4;
using namespace std::chrono;

void load_tles(const std::string& filename, std::vector<Tle>& tles) {
    std::ifstream file(filename);
    std::string line;
    std::string line1, line2;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line[0] == '1') {
            line1 = line.substr(0, 69);
        } else if (line[0] == '2') {
            line2 = line.substr(0, 69);
            try {
                tles.emplace_back("Benchmark", line1, line2);
            } catch (...) {}
        }
    }
}

int main(int argc, char* argv[]) {
    std::vector<Tle> tles;
    load_tles("SGP4-VER.TLE", tles);

    if (tles.empty()) {
        std::cerr << "No TLEs loaded!" << std::endl;
        return 1;
    }

    // Scale up the satellite count for Multi-Core test
    std::vector<Tle> scaled_tles;
    for(int i=0; i<300; i++) scaled_tles.insert(scaled_tles.end(), tles.begin(), tles.end());
    tles = scaled_tles;

    int iterations = 1000;
    if (argc > 1) iterations = std::stoi(argv[1]);

    std::cout << "Benchmarking " << tles.size() << " satellites and " << iterations << " iterations each..." << std::endl;

    // Scalar Benchmark
    {
        auto start = high_resolution_clock::now();
        unsigned long long total_propagations = 0;
        #pragma omp parallel for reduction(+:total_propagations)
        for (int t = 0; t < (int)tles.size(); ++t) {
            SGP4 sgp4(tles[t]);
            for (int i = 0; i < iterations; ++i) {
                try {
                    sgp4.FindPosition(static_cast<double>(i));
                    total_propagations++;
                } catch (...) { break; }
            }
        }
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(end - start).count();
        std::cout << "[Scalar] Time: " << duration/1e6 << "s, Avg: " << (double)duration/total_propagations << "us" << std::endl;
    }

    // Batch Benchmark
    {
        SGP4Batch batch(tles);
        std::vector<Eci> results;
        auto start = high_resolution_clock::now();
        unsigned long long total_propagations = 0;
        for (int i = 0; i < iterations; ++i) {
            batch.Propagate(static_cast<double>(i), results);
            for (const auto& eci : results) {
                if (eci.Position().x != 0.0) total_propagations++;
            }
        }
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(end - start).count();
        std::cout << "[Batch]  Time: " << duration/1e6 << "s, Avg: " << (double)duration/total_propagations << "us" << std::endl;
    }

    // Verification
    {
        SGP4Batch batch(tles);
        std::vector<Eci> results;
        batch.Propagate(10.0, results);
        int checks = 0;
        int failures = 0;
        for (size_t i = 0; i < tles.size(); ++i) {
            SGP4 sgp4(tles[i]);
            try {
                Eci expected = sgp4.FindPosition(10.0);
                OrbitalElements el(tles[i]);
                if (el.Period() < 225.0) {
                    checks++;
                    double diff = std::abs(expected.Position().x - results[i].Position().x);
                    if (diff > 1e-6) {
                        std::cout << "Mismatch for satellite " << i << ": diff=" << diff << std::endl;
                        failures++;
                    }
                }
            } catch(...) {}
        }
        std::cout << "[Verify] Accuracy check: " << checks << " objects verified, " << failures << " failures." << std::endl;
    }

    return 0;
}
