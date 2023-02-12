#pragma once

#include <fstream>
#include <vector>
#include <mutex>
#include <unordered_map>
#include "util/constants.h"
#include "util/basics.h"

namespace boltzi {
class BoundaryCurrentLog {
    private:
    std::mutex lock;
    std::vector<double> log;
    std::unordered_map<int, unsigned> origins;
    unsigned n_spawned = 0;
    #ifdef BOLTZI_BOUNDARY_DIRS
    std::vector<vec3> transition_dirs;
    #endif
    double t0;
    double dt;


    public:
    BoundaryCurrentLog(double t1, double t0, double dt) : log(static_cast<size_t>((t1-t0)/dt)), t0(t0), dt(dt) {}
    BoundaryCurrentLog(const BoundaryCurrentLog& other) {
        std::copy(other.log.begin(), other.log.end(), std::back_inserter(log));
        #ifdef BOLTZI_BOUNDARY_DIRS
        std::copy(other.transition_dirs.begin(), other.transition_dirs.end(), std::back_inserter(transition_dirs));
        #endif
        t0 = other.t0;
        dt = other.dt;
    }

    void transition(double t, double energy, vec3 dir, int source_idx, bool spawned = false) {
        std::lock_guard<std::mutex> guard(lock);
        log[static_cast<std::size_t>((t-t0)/dt)] += energy;
        #ifdef BOLTZI_BOUNDARY_DIRS
        transition_dirs.emplace_back(std::move(dir));
        #endif
        if(spawned)
            n_spawned += 1;
        else
            origins[source_idx] += 1;
    }
    void save(std::string name, double n_eff, double area) {
        std::lock_guard<std::mutex> guard(lock);
        std::ofstream file(name);
        for(unsigned i = 0; i<log.size();++i) {
            file<<i*dt+t0<<","<<log[i] * n_eff / area * util::HBAR / dt<<"\n";
        }
        file.close();

        #ifdef BOLTZI_BOUNDARY_DIRS
        std::ofstream file2(name+std::string(".dirs"));
        for(auto& entry : transition_dirs)
            file2<<entry[0]<<","<<entry[1]<<","<<entry[2]<<"\n";
        file2.close();
        #endif

        std::ofstream file3(name + std::string(".origins"));
        file3 << "spawn," << n_spawned << "\n";
        for(auto [source_idx, amount] : origins)
            file3 << source_idx << "," << amount << "\n";
        file3.close();
    }
};
}
