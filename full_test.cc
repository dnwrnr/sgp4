#include <SGP4.h>
#include <Tle.h>
#include <Vector.h>
#include <DateTime.h>
#include <SolarPosition.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <map>
#include "rf_model_full.h"

struct SolarData {
    std::map<long long, std::vector<double>> data;
};

SolarData LoadSolarFull(const std::string& path) {
    SolarData data;
    std::ifstream file(path);
    std::string line;
    std::getline(file, line); // header
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string val;
        std::vector<std::string> row;
        while(std::getline(ss, val, ',')) row.push_back(val);
        if (row.size() >= 28) {
            try {
                long long key = std::stoll(row[0].substr(0,4) + row[0].substr(5,2) + row[0].substr(8,2));
                data.data[key] = {std::stod(row[24]), std::stod(row[27]), std::stod(row[20])};
            } catch (...) { continue; }
        }
    }
    return data;
}

int main(int argc, char* argv[]) {
    if (argc < 3) return 1;
    SolarData solar = LoadSolarFull(argv[2]);
    std::ifstream hfile(argv[1]);
    std::string l1, l2;
    std::vector<libsgp4::Tle> history;
    while(std::getline(hfile, l1) && std::getline(hfile, l2)) {
        if (l1.length() < 69 || l2.length() < 69) continue;
        history.push_back(libsgp4::Tle(l1.substr(0, 69), l2.substr(0, 69)));
    }
    if (history.empty()) return 0;

    libsgp4::Tle startTle = history[0];
    libsgp4::SGP4 propagator(startTle);
    double initial_bstar = startTle.BStar();
    double k_factor = 0.8; 

    std::cout << "Days,Residual_km,Bias_rad" << std::endl;

    for (size_t i = 1; i < history.size(); ++i) {
        double minutes = (history[i].Epoch() - startTle.Epoch()).TotalMinutes();
        double days_elapsed = minutes / 1440.0;
        
        double sum_delta_b = 0;
        double sum_sum_delta_b = 0;
        for (int d = 0; d <= (int)days_elapsed; ++d) {
            std::string dstr = startTle.Epoch().AddDays(d-1).ToString();
            long long key = std::stoll(dstr.substr(0,4) + dstr.substr(5,2) + dstr.substr(8,2));
            double f107 = (solar.data.count(key)) ? solar.data[key][0] : 70;
            double b_pred = PredictBStarFull(f107, 70, 5, 30.0);
            sum_delta_b += (b_pred - initial_bstar);
            sum_sum_delta_b += sum_delta_b;
        }
        
        propagator.SetAlongTrackBias(k_factor * sum_sum_delta_b);
        
        libsgp4::Eci predicted = propagator.FindPosition(minutes);
        libsgp4::Eci actual = libsgp4::SGP4(history[i]).FindPosition(0.0);
        
        std::cout << days_elapsed << "," << (predicted.Position() - actual.Position()).Magnitude() << "," << k_factor * sum_sum_delta_b << std::endl;
    }
    return 0;
}
