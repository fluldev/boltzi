#pragma once

#include <string>
#include<vector>
#include<queue>
#include <mutex>
#include <fstream>

#include "boltzi_exception.h"
#include "boltzi_types.h"

namespace boltzi {
    class PhononData {
        public:
        vec3 k;
        vec3 r;
        unsigned material_idx;
        unsigned branch_idx;
        short sign;
        bool active;
        double last_event_path = 0;
        double last_event_time = 0;
        double init_time;
        int source_idx;  // -1 = inital volume phonons, -2 = spatial phonons, >=0 = boundary_idx

        PhononData() = delete;
        PhononData(vec3 k, vec3 r, unsigned material_idx, unsigned branch_idx, short sign, double init_time, int source_idx = -1)
            : k(k), r(r), material_idx(material_idx), branch_idx(branch_idx), sign(sign), init_time(init_time), source_idx(source_idx) {
                active = true;
                last_event_time = 0;
                last_event_path = 0;
            }
    };

    class PhononList {
        private:
            std::queue<unsigned> free_spots{};
            std::mutex free_spots_lock{};

        public:
            using phonon_iterator = std::vector<PhononData>::iterator;
            std::vector<PhononData> phonons{};
            void apply(void (*fn)(phonon_iterator phonon)) {
                // simple implementation for now should be run in parallel
                for(auto it = phonons.begin(); it!=phonons.end(); ++it) {
                    if(it->active) fn(it);
                }
            }

            void clear() {
                free_spots = std::queue<unsigned>{};
                phonons.clear();
            }
            
            void add(const PhononData& phonon) {
                // this should not be used in parallel code
                if(!free_spots.empty()) {
                    phonons[free_spots.front()] = phonon;
                    free_spots.pop(); 
                }
                else {
                    std::lock_guard<std::mutex> cur_lock(free_spots_lock);
                    phonons.push_back(phonon);
                }
            }
            
            void remove(phonon_iterator iterator) {
                if(!iterator->active) return;  // already removed
                std::lock_guard<std::mutex> cur_lock(free_spots_lock);
                iterator->active = false;
                free_spots.push(iterator-phonons.cbegin());
            }

            int n_active() const {
                return phonons.size() - free_spots.size();
            }
            void save(std::string name) const {
                std::ofstream file(name);
                file<<"simulation_state"<<std::endl;
                for(auto& phonon : phonons)
                    if(phonon.active) {
                        file<<phonon.k[0]<<","<<phonon.k[1]<<","<<phonon.k[2]<<","; 
                        file<<phonon.r[0]<<","<<phonon.r[1]<<","<<phonon.r[2]<<","; 
                        file<<phonon.material_idx<<","<<phonon.branch_idx<<"\n";
                        file<<phonon.sign<<","<<phonon.active<<"\n";                        
                        file<<phonon.last_event_path<<","<<phonon.last_event_time<<","<<phonon.init_time<<"\n";
                    }
                file.close();
            }

            void load(std::string name) {
                std::ifstream file(name);
                if(!file.is_open()) throw BoltziException("Can't open file.");
                clear();
                std::string line;
                file>>line;  // ignore first line
                file>>line;  // ignore first line
                size_t cursor=0;

                auto get_one = [&]()->std::string {
                    size_t pos = line.find(",", cursor);
                    if(pos == std::string::npos) {
                        std::string res = line.substr(cursor, line.length()-cursor);
                        cursor = 0;
                        file>>line;
                        return res;
                    }
                    unsigned old_cursor = cursor;
                    cursor = pos+1;
                    return line.substr(old_cursor, pos-old_cursor);
                };

                try {
                    while(!file.eof()) {
                        PhononData phonon{vec3{0,0,0}, vec3{0,0,0}, 0, 0, 0, 0};
                        phonon.k[0] = std::stod(get_one());
                        phonon.k[1] = std::stod(get_one());
                        phonon.k[2] = std::stod(get_one());

                        phonon.r[0] = std::stod(get_one());
                        phonon.r[1] = std::stod(get_one());
                        phonon.r[2] = std::stod(get_one());

                        phonon.material_idx = std::stoi(get_one());
                        phonon.branch_idx = std::stoi(get_one());
                        phonon.sign = std::stoi(get_one());
                        phonon.active = std::stoi(get_one());
                        phonon.last_event_path = std::stod(get_one());
                        phonon.last_event_time = std::stod(get_one());
                        phonon.init_time = std::stod(get_one());

                        phonons.push_back(phonon);
                        if(!phonon.active)
                            free_spots.push(phonons.size()-1);
                    }
                }
                catch(...) {throw BoltziException("Save file is corrupted.");}
                file.close();
        }
    };
}
