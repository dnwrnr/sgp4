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
#include "DateTime.h"

using namespace libsgp4;
using namespace std;

struct SatState {
    string id;
    string l1, l2;
};

int main(int argc, char* argv[]) {
    vector<SatState> test_sats = {
        {"ISS (LEO)", 
         "1 25544U 98067A   24001.01267188  .00016541  00000-0  29758-3 0  5694", 
         "2 25544  51.6422  68.6294 0003347 343.4617  78.0593 15.49961425432470"},
        {"GPS (MEO)",
         "1 43637U 18078A   24001.44203794  .00000031  00000-0  00000+0 0  9990",
         "2 43637  55.3371  53.5134 0003463 269.9676  90.0076  2.00557457 37777"},
        {"GOES 16 (GEO)",
         "1 41866U 16071A   24001.40149023  .00000084  00000-0  00000+0 0  9994",
         "2 41866   0.0242 101.4390 0001043 273.7380 344.6041  1.00273295 25958"}
    };

    cout << left << setw(15) << "Object" << setw(15) << "Mode" << setw(15) << "Time (min)" 
         << setw(20) << "Standard Pos (km)" << setw(15) << "Error (m)" << endl;
    cout << string(85, '-') << endl;

    for (const auto& s : test_sats) {
        Tle tle(s.l1, s.l2);
        SGP4 model(tle);
        
        vector<Tle> tles = {tle};
        SGP4Batch batch(tles);
        vector<Eci> results_std, results_fast;

        double test_times[] = {0.0, 1000.0, 10000.0};
        for (double tsince : test_times) {
            batch.Propagate(tsince, results_std, SGP4Batch::MathMode::Standard);
            batch.Propagate(tsince, results_fast, SGP4Batch::MathMode::Minimax);
            
            // Standard Reference (Single Object)
            Vector std_pos = model.FindPosition(tsince).Position();
            
            // Batch Standard (SIMD Vector Math)
            Vector batch_std = results_std[0].Position();
            
            // Batch Fast (Custom Minimax)
            Vector batch_fast = results_fast[0].Position();

            double err_std = (std_pos - batch_std).Magnitude();
            double err_fast = (std_pos - batch_fast).Magnitude();
            
            cout << left << setw(15) << s.id << setw(15) << "Standard" << setw(15) << tsince 
                 << fixed << setprecision(2) << setw(20) << std_pos.Magnitude() 
                 << scientific << setprecision(4) << err_std * 1000.0 << endl;
            cout << left << setw(15) << "" << setw(15) << "Minimax" << setw(15) << tsince 
                 << fixed << setprecision(2) << setw(20) << batch_fast.Magnitude() 
                 << scientific << setprecision(4) << err_fast * 1000.0 << endl;
        }
        cout << endl;
    }

    // Performance Test
    vector<Tle> many_tles;
    for(int i=0; i<10000; i++) many_tles.push_back(Tle(test_sats[0].l1, test_sats[0].l2));
    SGP4Batch perf_batch(many_tles);
    vector<Eci> res;
    
    // Warmup
    for(int t=0; t<10; t++) perf_batch.Propagate(t, res, SGP4Batch::MathMode::Standard);

    auto start = chrono::high_resolution_clock::now();
    for(int t=0; t<100; t++) perf_batch.Propagate(t, res, SGP4Batch::MathMode::Standard);
    auto end = chrono::high_resolution_clock::now();
    double time_std = chrono::duration<double>(end-start).count();

    start = chrono::high_resolution_clock::now();
    for(int t=0; t<100; t++) perf_batch.Propagate(t, res, SGP4Batch::MathMode::Minimax);
    end = chrono::high_resolution_clock::now();
    double time_fast = chrono::duration<double>(end-start).count();

    cout << "Throughput Comparison (10k objects, 100 iterations):" << endl;
    cout << "Standard (libmvec): " << (1000000.0 / time_std / 1e6) << " million sats/sec" << endl;
    cout << "Minimax (Custom):   " << (1000000.0 / time_fast / 1e6) << " million sats/sec" << endl;

    return 0;
}
