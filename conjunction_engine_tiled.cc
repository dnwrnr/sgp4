#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <algorithm>
#include <immintrin.h>
#include "Tle.h"
#include "SGP4.h"
#include "SGP4Batch.h"
#include "Globals.h"
#include "Eci.h"

using namespace libsgp4;
using namespace std;

struct RiskData {
    double volatility = 0.0;
    int maneuver_count = 0;
};

struct Conjunction {
    int id1;
    int id2;
    double tsince;
    float miss_distance;
    float prox_index;
    float volatility_buf;
};

map<int, RiskData> load_global_risk_data() {
    map<int, RiskData> risk;
    ifstream tfile("tumble_residual_results.csv");
    string line;
    if (tfile.is_open()) {
        getline(tfile, line);
        while (getline(tfile, line)) {
            stringstream ss(line);
            string sid, vol;
            getline(ss, sid, ',');
            getline(ss, vol, ',');
            try { if (!sid.empty()) risk[stoi(sid)].volatility = stod(vol); } catch(...) {}
        }
        tfile.close();
    }
    ifstream mfile("maneuvers_report.csv");
    if (mfile.is_open()) {
        getline(mfile, line);
        while (getline(mfile, line)) {
            stringstream ss(line);
            string sid, mcount;
            getline(ss, sid, ',');
            getline(ss, mcount, ',');
            try { if (!sid.empty()) risk[stoi(sid)].maneuver_count = stoi(mcount); } catch(...) {}
        }
        mfile.close();
    }
    return risk;
}

struct SatelliteInfo {
    int id;
    float risk_vol;
    int risk_man;
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <tle_file> [limit]" << endl;
        return 1;
    }

    map<int, RiskData> global_risk = load_global_risk_data();
    vector<Tle> tles;
    vector<SatelliteInfo> sat_info;

    cout << "Loading Satellites from " << argv[1] << "..." << endl;
    ifstream tfile(argv[1]);
    string line;
    int limit = (argc > 2) ? stoi(argv[2]) : 10000;
    
    map<int, pair<string, string>> latest_tle_strings;
    string l1, l2;
    while (getline(tfile, line)) {
        if (line.empty()) continue;
        if (line.back() == '\\' || line.back() == '\r') line.pop_back();
        
        if (line[0] == '1') {
            l1 = line;
        } else if (line[0] == '2') {
            l2 = line;
            if (l1.empty()) continue;
            try {
                int id = stoi(l1.substr(2, 5));
                if (latest_tle_strings.find(id) == latest_tle_strings.end()) {
                    if (latest_tle_strings.size() < (size_t)limit) {
                         latest_tle_strings[id] = {l1, l2};
                    }
                }
            } catch(...) {}
            l1.clear();
        }
    }
    tfile.close();

    for (auto const& [id, lines] : latest_tle_strings) {
        try {
            string s1 = lines.first.substr(0, 69);
            string s2 = lines.second.substr(0, 69);
            Tle tle(s1, s2);
            tles.push_back(tle);
            sat_info.push_back({id, (float)global_risk[id].volatility, global_risk[id].maneuver_count});
        } catch(...) {}
    }

    cout << "Initializing Batch Propagator for " << tles.size() << " satellites..." << endl;
    if (tles.empty()) {
        cerr << "Error: No satellites loaded." << endl;
        return 1;
    }
    SGP4Batch batch(tles);
    int n = tles.size();

    alignas(64) vector<double> pos_x(n), pos_y(n), pos_z(n);
    vector<Eci> results;

    cout << "Performing Tiled Screening (1-hour window, 10km filter)..." << endl;
    vector<Conjunction> events;
    const double sigma_sys = 1.0;
    const double sq_threshold = 100.0;
    const int TILE_SIZE = 256; 

    for (int t = 0; t < 60; ++t) {
        batch.Propagate((double)t, results);
        
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            Vector p = results[i].Position();
            pos_x[i] = p.x; pos_y[i] = p.y; pos_z[i] = p.z;
        }

        #pragma omp parallel
        {
            vector<Conjunction> local_events;
            local_events.reserve(100);
            #pragma omp for schedule(dynamic)
            for (int ii = 0; ii < n; ii += TILE_SIZE) {
                for (int jj = ii; jj < n; jj += TILE_SIZE) {
                    int i_end = min(ii + TILE_SIZE, n);
                    int j_end = min(jj + TILE_SIZE, n);

                    for (int i = ii; i < i_end; i++) {
                        double xi = pos_x[i], yi = pos_y[i], zi = pos_z[i];
                        int start_j = (ii == jj) ? (i + 1) : jj;

                        int j = start_j;
                        __m512d v_xi = _mm512_set1_pd(xi);
                        __m512d v_yi = _mm512_set1_pd(yi);
                        __m512d v_zi = _mm512_set1_pd(zi);
                        __m512d v_thresh = _mm512_set1_pd(sq_threshold);

                        for (; j <= j_end - 8; j += 8) {
                            __m512d v_xj = _mm512_loadu_pd(&pos_x[j]);
                            __m512d v_yj = _mm512_loadu_pd(&pos_y[j]);
                            __m512d v_zj = _mm512_loadu_pd(&pos_z[j]);

                            __m512d dx = _mm512_sub_pd(v_xi, v_xj);
                            __m512d dy = _mm512_sub_pd(v_yi, v_yj);
                            __m512d dz = _mm512_sub_pd(v_zi, v_zj);

                            __m512d d2 = _mm512_mul_pd(dx, dx);
                            d2 = _mm512_fmadd_pd(dy, dy, d2);
                            d2 = _mm512_fmadd_pd(dz, dz, d2);

                            __mmask8 mask = _mm512_cmp_pd_mask(d2, v_thresh, _CMP_LT_OQ);
                            if (mask != 0) {
                                alignas(64) double dists_sq[8];
                                _mm512_store_pd(dists_sq, d2);
                                for (int k = 0; k < 8; k++) {
                                    if (mask & (1 << k)) {
                                        float d_sq = (float)dists_sq[k];
                                        Conjunction c;
                                        c.id1 = sat_info[i].id;
                                        c.id2 = sat_info[j+k].id;
                                        c.tsince = t;
                                        c.miss_distance = sqrtf(d_sq);
                                        c.prox_index = expf(-(d_sq) / (float)(2.0 * sigma_sys * sigma_sys));
                                        float v1 = 1.0f + (sat_info[i].risk_man > 0 ? 1.0f : 0.0f) + sat_info[i].risk_vol;
                                        float v2 = 1.0f + (sat_info[j+k].risk_man > 0 ? 1.0f : 0.0f) + sat_info[j+k].risk_vol;
                                        c.volatility_buf = max(v1, v2);
                                        local_events.push_back(c);
                                    }
                                }
                            }
                        }
                        for (; j < j_end; j++) {
                            double dx = xi - pos_x[j], dy = yi - pos_y[j], dz = zi - pos_z[j];
                            double d2 = dx*dx + dy*dy + dz*dz;
                            if (d2 < sq_threshold) {
                                float d_sq = (float)d2;
                                Conjunction c;
                                c.id1 = sat_info[i].id;
                                c.id2 = sat_info[j].id;
                                c.tsince = t;
                                c.miss_distance = sqrtf(d_sq);
                                c.prox_index = expf(-(d_sq) / (float)(2.0 * sigma_sys * sigma_sys));
                                float v1 = 1.0f + (sat_info[i].risk_man > 0 ? 1.0f : 0.0f) + sat_info[i].risk_vol;
                                float v2 = 1.0f + (sat_info[j].risk_man > 0 ? 1.0f : 0.0f) + sat_info[j].risk_vol;
                                c.volatility_buf = max(v1, v2);
                                local_events.push_back(c);
                            }
                        }
                    }
                }
            }
            #pragma omp critical
            {
                if (events.size() < 1000000) {
                    events.insert(events.end(), local_events.begin(), local_events.end());
                }
            }
         local_events.clear();
        }
    }

    sort(events.begin(), events.end(), [](const Conjunction& a, const Conjunction& b) {
        return a.prox_index > b.prox_index;
    });

    cout << "Writing Dual-Vector Risk Report: dual_vector_report.csv" << endl;
    ofstream outfile("dual_vector_report.csv");
    outfile << "ID1,ID2,T_Since,Miss_Distance_km,Proximity_Index,Volatility_Buffer,Quadrant" << endl;
    int count = 0;
    for (const auto& e : events) {
        if (++count > 100000) break;
        string quadrant;
        if (e.prox_index > 0.5f && e.volatility_buf < 2.0f) quadrant = "Deterministic";
        else if (e.prox_index < 0.5f && e.volatility_buf > 2.0f) quadrant = "Stochastic";
        else if (e.prox_index > 0.5f && e.volatility_buf > 2.0f) quadrant = "Critical";
        else quadrant = "Background";

        outfile << e.id1 << "," << e.id2 << "," << e.tsince << "," << fixed << setprecision(3)
                << e.miss_distance << "," << e.prox_index << "," << e.volatility_buf << "," << quadrant << endl;
    }

    cout << "Analysis Complete. Found " << events.size() << " proximity events." << endl;
    return 0;
}
