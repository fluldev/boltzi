//#include <bits/chrono.h>
#include <cmath>
#include <limits>
#include <utility>
#include <memory>
#include <unordered_map>
#include <armadillo>
#include <atomic>
#include <chrono>

#include <iostream>

#include "headers/material.h"
#include "headers/phonon.h"
#include "headers/sample.h"
#include "headers/boltzi_config.h"
#include "headers/boltzi_types.h"
#include "headers/util/pool.h"
#include "headers/util/random.h"
#include "headers/silencer.h"
#include "headers/boundary_current_log.h"

using boltzi::config::Silencer;

void boltzi::Sample::calc_bin_temps() {
    struct BranchInfo {
        // all members do not get written to simultaneously by multiple threads since 
        // all threads on this struct work on diffrent cells
        volatile unsigned n_phonons{0};  // amount of phonons of this branch
        volatile double omega_sum{0.0};  // plain omega sum
        volatile double* scatter_sums; 
        BranchInfo() : n_phonons{0}, omega_sum{0.0}, scatter_sums{nullptr} {}
    };


    // Indexed by cell ijk, material_idx, branch_idx (not all cells contain all materials 
    // and not all cells contain phonons from all branches of a specifc material)
    util::Grid<3, 
        std::unordered_map<unsigned int, std::unordered_map<unsigned int, BranchInfo>>> 
                omega_sums{std::get<0>(n_bins), std::get<1>(n_bins), std::get<2>(n_bins)};

    //TODO: reallocating all the time is inefficient better just reuse omegasums and set everything to zero before (maybe in multiple threads)
    #ifdef BOLTZI_MEASURE_TIMINGS
    std::cout<<"cell_temp timings: ";
    #endif
    auto t0 = std::chrono::high_resolution_clock::now();
    // calculate omega sums
    for(unsigned phonon_idx = 0; phonon_idx < phonons.phonons.size(); ++phonon_idx) {
        auto& phonon = phonons.phonons[phonon_idx];
        if(!phonon.active) continue;
        if(!this->geometry.is_inside_box(phonon.r)) {
            this->phonons.remove(phonon_idx + phonons.phonons.begin());  // this doesn't change the phonon vector and therefore doesn't invalidate any iterator
            continue;
        }
        auto[i,j,k] = this->pos_to_cell(phonon.r);
        // TODO This access could add Elements to the unordered_map and is therefore not thread safe
        // access with at would be thread-save since the elements are atomics and at doesn't add new elements to the unordered_map
        auto& cur = omega_sums.at({i,j,k})[phonon.material_idx][phonon.branch_idx];
        //if(cur.scatter_sums == nullptr) {
        //    unsigned n_three_prc = materials[phonon.material_idx].branches[phonon.branch_idx]->three_phonon_processes.size();
        //    if(n_three_prc>0)
        //        cur.scatter_sums = new double[n_three_prc];
        //}
        cur.n_phonons += 1;
        cur.omega_sum += omega(phonon) * phonon.sign;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    #ifdef BOLTZI_MEASURE_TIMINGS
    std::cout<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0)<<", ";
    #endif

    // calculate temperatures
    // dangerous code maybe since omega_sums elements get modified by multiple threads (however not the same by multiple threads at once)
    t0 = std::chrono::high_resolution_clock::now();
    util::pool.add_container_job(
        [this] (auto cell_it) {
            double weight_sum=0;
            double temp_sum=0;
            auto [i,j,k] = (*cell_it).first;
            // loops go only over occupied branches (phonons of this branch exist in cell ijk)
            for(auto material_it = (*cell_it).second.begin(); material_it != (*cell_it).second.end(); ++material_it) {
                auto material_idx = material_it->first;
                auto& material = material_it->second;
                for(auto branch_it = material.begin(); branch_it!=material.end();++branch_it) {
                    auto branch_idx = branch_it->first;
                    auto& branch_info = branch_it->second;
                    unsigned n_three_prc = materials[material_idx].branches[branch_idx]->three_phonon_processes.size();
                    if(n_three_prc>0)
                        branch_info.scatter_sums = new double[n_three_prc];
                    for(unsigned i = 0; i<n_three_prc;++i)
                        branch_info.scatter_sums[i] = 0.0;
                    //double cur_weight = materials[material_idx].branches[branch_idx]->branch_weight() * branch_info.n_phonons;  // TODO: is this the better choice?
                    double cur_weight = branch_info.n_phonons;
                    if (branch_info.n_phonons >= config::MINIMUM_TEMPERATURE_PHONONS) {
                        double cur_temp = this->materials[material_idx].branches[branch_idx]->omegasum_to_temp(
                            branch_info.omega_sum / this->cell_volumes.at({i,j,k})[material_idx] * this->n_eff
                            + this->materials[material_idx].branches[branch_idx]->calc_omegasum(this->deviational_temp(this->cell_to_pos(i, j, k)))
                        );
                        if(!std::isnan(cur_temp))
                            temp_sum += cur_temp * cur_weight;
                        else {
                            cur_weight = config::MINIMUM_TEMPERATURE_PHONONS;
                            temp_sum += this->deviational_temp(this->cell_to_pos(i, j, k)) * cur_weight;
                        }
                    }
                    else {
                        cur_weight = config::MINIMUM_TEMPERATURE_PHONONS;
                        temp_sum+=this->deviational_temp(this->cell_to_pos(i, j, k)) * cur_weight;
                    }
                    weight_sum += cur_weight;
                }
            }
            if(weight_sum>0)
                this->cell_temps.at({i,j,k}) = temp_sum / weight_sum;
            else
                this->cell_temps.at({i,j,k}) = this->deviational_temp(this->cell_to_pos(i, j, k));
        },
        omega_sums
    );
    util::pool.exec();
    t1 = std::chrono::high_resolution_clock::now();
    #ifdef BOLTZI_MEASURE_TIMINGS
    std::cout<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0)<<", ";
    #endif
    
    // calculate scatter sums
    t0 = std::chrono::high_resolution_clock::now();
    for(unsigned phonon_idx = 0; phonon_idx < phonons.phonons.size(); ++phonon_idx) {
        auto& phonon = phonons.phonons[phonon_idx];
        if(!phonon.active) continue;
        auto[i,j,k] = pos_to_cell(phonon.r);
        auto& cur = omega_sums.at({i,j,k})[phonon.material_idx][phonon.branch_idx];
        unsigned n_three_prc = materials[phonon.material_idx].branches[phonon.branch_idx]->three_phonon_processes.size();
        for(unsigned prc_idx = 0; prc_idx<n_three_prc;++prc_idx)
            cur.scatter_sums[prc_idx] += omega(phonon) * phonon.sign\
            * materials[phonon.material_idx].branches[phonon.branch_idx]->three_phonon_processes[prc_idx]->scatter_rate(phonon, cell_temps.at({i,j,k}));
    }
    t1 = std::chrono::high_resolution_clock::now();
    #ifdef BOLTZI_MEASURE_TIMINGS
    std::cout<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0)<<", ";
    #endif

    // calculate pseudo temperatures
    t0 = std::chrono::high_resolution_clock::now();
    util::pool.add_container_job(
        [this] (const auto cell_it) {
            auto [i,j,k] = (*cell_it).first;
            for(auto material_it = (*cell_it).second.begin(); material_it != (*cell_it).second.end(); ++material_it) {
                auto material_idx = material_it->first;
                auto& material = material_it->second;
                for(auto branch_it = material.begin(); branch_it!=material.end();++branch_it) {
                    auto branch_idx = branch_it->first;
                    auto& branch_info = branch_it->second;
                    if(branch_info.n_phonons == 0) continue;  // will not be used by any phonons -> does not need to be calculated
                                                              // there was also no branch_info.scatter_sums memory allocated

                    unsigned n_three_prc = materials[material_idx].branches[branch_idx]->three_phonon_processes.size();
                    if (branch_info.n_phonons >= config::MINIMUM_TEMPERATURE_PHONONS) {
                        for(unsigned prc_idx = 0; prc_idx<n_three_prc;++prc_idx) {
                            ThreePhononProcess& prc = *(materials[material_idx].branches[branch_idx]->three_phonon_processes[prc_idx]);

                            double cur_pt = prc.pseudo_temp(
                                branch_info.scatter_sums[prc_idx] / cell_volumes.at({i,j,k})[material_idx] * n_eff\
                                + prc.scattered_omegasum(cell_temps.at({i,j,k}), deviational_temp(cell_to_pos(i,j,k)))
                                , cell_temps.at({i,j,k})  // might be better to use the temperature 
                                                          // of this specifc branch only but this would require even more memory
                                                          // the temperature of a specific branch doesn't make sense however
                            );
                            if(!std::isnan(cur_pt))
                                cell_pseudo_temps.at({i,j,k})[material_idx][branch_idx][prc_idx] = cur_pt;
                            else
                                cell_pseudo_temps.at({i,j,k})[material_idx][branch_idx][prc_idx] = cell_temps.at({i,j,k}); 
                        }
                    }
                    else
                        for(unsigned prc_idx = 0; prc_idx<n_three_prc; ++prc_idx) 
                            cell_pseudo_temps.at({i,j,k})[material_idx][branch_idx][prc_idx] = cell_temps.at({i,j,k});
                    if(n_three_prc>0)
                        delete [] branch_info.scatter_sums;
                }
            }
        },
        omega_sums
    );
    util::pool.exec();
    t1 = std::chrono::high_resolution_clock::now();
    #ifdef BOLTZI_MEASURE_TIMINGS
    std::cout<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0)<<", ";
    std::cout<<std::endl;
    #endif
}


boltzi::Sample::Sample(
    double n_eff,
    std::tuple<unsigned, unsigned, unsigned> n_bins,
    const SampleGeometry& geometry,
    const std::vector<Material>& materials,
    std::function<double(vec3 r)> deviational_temp,
    std::function<vec3(vec3 r)> deviational_temp_gradient 
) : n_bins(n_bins),
    cell_size{
        (geometry.box.second[0]-geometry.box.first[0]) / std::get<0>(n_bins),
        (geometry.box.second[1]-geometry.box.first[1]) / std::get<1>(n_bins),
        (geometry.box.second[2]-geometry.box.first[2]) / std::get<2>(n_bins)
    },
    cell_temps{std::get<0>(n_bins), std::get<1>(n_bins), std::get<2>(n_bins)}, 
    cell_pseudo_temps{std::get<0>(n_bins), std::get<1>(n_bins), std::get<2>(n_bins)}, 
    materials(materials), 
    n_eff(n_eff),
    cell_volumes{std::get<0>(n_bins), std::get<1>(n_bins), std::get<2>(n_bins)}, 
    deviational_temp(deviational_temp),
    deviational_temp_gradient(deviational_temp_gradient),
    geometry(geometry)
{
    // initialize pseudo_temp arrays
    for(unsigned i = 0; i<std::get<0>(n_bins); ++i) {
        for(unsigned j = 0; j<std::get<1>(n_bins); ++j) {
            for(unsigned k = 0; k<std::get<2>(n_bins); ++k) {
                unsigned material_idx = 0; 
                for(const auto& material : materials) {
                    unsigned branch_idx = 0;
                    for(const auto& branch : material.branches) {
                        cell_pseudo_temps.at({i,j,k})[material_idx][branch_idx] = new double [branch->three_phonon_processes.size()];
                        ++branch_idx;
                    }
                    ++material_idx;
                }
            }}}

    }

boltzi::Sample::~Sample() {
    for(unsigned i = 0; i<std::get<0>(n_bins); ++i) {
        for(unsigned j = 0; j<std::get<1>(n_bins); ++j) {
            for(unsigned k = 0; k<std::get<2>(n_bins); ++k) {
                unsigned material_idx = 0; 
                for(const auto& material : materials) {
                    for(unsigned branch_idx=0;branch_idx<material.branches.size(); ++branch_idx) {
                        delete [] cell_pseudo_temps.at({i,j,k})[material_idx][branch_idx];
                    }
                    ++material_idx;
                }
            }}}

}


void boltzi::Sample::init_cells() {
    vec3 cell_unit =    (geometry.box.second-geometry.box.first) 
                        / vec3{
                            static_cast<double>(std::get<0>(n_bins)), 
                            static_cast<double>(std::get<1>(n_bins)), 
                            static_cast<double>(std::get<2>(n_bins))
                        };
    double cell_volume_unit = arma::prod(cell_unit) / config::CELL_OVERLAP_SAMPLES;
    util::pool.add_container_job(
        [this,cell_unit,cell_volume_unit](decltype(cell_volumes)::iterator it) {
            auto [i,j,k] = (*it).first;
            for(unsigned material_idx = 0;material_idx<materials.size();++material_idx) 
                (*it).second[material_idx] = 0;

            for(unsigned n = 0; n<config::CELL_OVERLAP_SAMPLES;++n) {
                vec3 r = util::random.random_pos_cell(
                    geometry.box.first 
                    + cell_unit % vec3{
                        static_cast<double>(i),
                        static_cast<double>(j),
                        static_cast<double>(k)
                    }, cell_unit);
                if(geometry.is_inside(r))
                    (*it).second[material_index(r)] += cell_volume_unit;               
            }
        }
    , cell_volumes);
    util::pool.exec();
}

unsigned boltzi::Sample::spawn_phonons(std::function<double(vec3)> init_temp, double init_time, bool use_turnover) {
    // calculate amount of required phonons using the cell volumes
    // must be a double in case of <1 phonons spawns for small cells, small dt or large n_eff
    std::unordered_map<unsigned, std::unordered_map<unsigned, double>> n_phonons;
    unsigned n_spawned=0;
    for(unsigned material_idx = 0; material_idx<materials.size();++material_idx)
        for(unsigned branch_idx = 0; branch_idx<materials[material_idx].branches.size(); ++branch_idx) {
            if(use_turnover)
                n_phonons[material_idx][branch_idx] = volume_turnover[material_idx][branch_idx];
            for(auto it = cell_volumes.begin(); it != cell_volumes.end(); ++it) {
                auto [i,j,k] = (*it).first;
                vec3 cellpos = cell_to_pos(i,j,k);
                n_phonons[material_idx][branch_idx] += materials[material_idx].branches[branch_idx]->phonon_density(
                        init_temp(cellpos), deviational_temp(cellpos)
                ) / n_eff * (*it).second[material_idx];
            }
            if(use_turnover)
                volume_turnover[material_idx][branch_idx] = util::modf(n_phonons[material_idx][branch_idx]);
        }
   

    unsigned material_idx = 0;
    for(const auto& material : materials) {
        unsigned branch_idx = 0;
        for(const auto& branch : material.branches) { 
            n_spawned+=n_phonons[material_idx][branch_idx];
            auto positions = geometry.random_positions(static_cast<unsigned>(n_phonons[material_idx][branch_idx]), 
                [this, &init_temp, &branch](vec3 pos)->double{
                    return branch->phonon_density(init_temp(pos), this->deviational_temp(pos));
                });

            for(auto pos : positions) {
                double temp = init_temp(pos);
                double dev_temp = deviational_temp(pos);
                auto k = branch->draw_k(temp, dev_temp);
                phonons.add(PhononData{
                    k,
                    pos,
                    material_idx,
                    branch_idx,
                    static_cast<short>((temp>=dev_temp)*2-1),
                    init_time,
                    -1
                });
            }
            ++branch_idx;
        }
        ++material_idx;
    }
    return n_spawned;
}

unsigned boltzi::Sample::spawn_surface_phonons(double dt, double init_time, bool use_turnover, std::unordered_map<unsigned, BoundaryCurrentLog>* collision_log) {
    unsigned n_spawned=0;
    unsigned boundary_idx = 0;
    for(const auto& boundary : geometry.boundaries) {
        if(!boundary->is_generator()) {++boundary_idx; continue;}
        double n_phonons = boundary->surface_integral([this, boundary] (vec3 pos) -> double {
            double sum = 0;
            for(const auto& branch : this->associated_material(pos).branches)
                sum +=  branch->phonon_area_flux(boundary->surface_temperature(pos), this->deviational_temp(pos));
            return sum;
        }) / n_eff * dt;
        if(use_turnover) {
            n_phonons += area_turnover[boundary_idx];
            area_turnover[boundary_idx] = util::modf(n_phonons);
        }
        n_spawned+=n_phonons;

        std::vector<std::tuple<vec3,vec3,vec3,vec3>> pos_and_surf{boundary->random_positions(
            static_cast<unsigned>(n_phonons),
            [this, boundary] (vec3 pos) -> double {
            double sum = 0;
            for(const auto& branch : this->associated_material(pos).branches)
                sum += branch->phonon_area_flux(boundary->surface_temperature(pos), this->deviational_temp(pos));
            return sum;
        })};

        for(const auto& [pos, normal, e1, e2] : pos_and_surf) {
            //TODO: doing this for each position might be a performance impact
            std::vector<double> branch_weights;
            double temp = boundary->surface_temperature(pos);
            double dev_temp = deviational_temp(pos);
            
            unsigned material_idx = material_index(pos);

            const auto& material = materials[material_idx];               
            for(const auto& branch : material.branches)
                branch_weights.push_back(branch->phonon_density(temp, dev_temp) * branch->avg_velocity());

            unsigned branch_idx = util::random.draw_index_weighted(branch_weights);
            auto& branch = material.branches[branch_idx];
            auto k = branch->draw_k_area(temp, dev_temp, normal, e1, e2);
            PhononData cur_phonon{
                k,
                pos,
                material_idx,
                branch_idx,
                static_cast<short>((temp>=dev_temp)*2-1),
                init_time + util::random.draw_real()*dt,
                static_cast<int>(boundary_idx)
            };
            if(collision_log)
                collision_log->at(boundary_idx).transition(
                    cur_phonon.init_time, 
                    -omega(cur_phonon) * cur_phonon.sign, 
                    cur_phonon.k * cur_phonon.sign,
                    cur_phonon.source_idx,
                    true
                );
            phonons.add(std::move(cur_phonon));
        }
        boundary_idx+=1;
    }
    return n_spawned;
}


unsigned boltzi::Sample::spawn_spatial_phonons(double dt, double init_time, bool use_turnover) {
    //TODO: Die Rate muss nur einmalig berechnet werden, da sich die deviational temp ja nicht ändert!
    // damit fliegt auch die config option raus
    // calculate amount of required phonons using the cell volumes
    // Selbiges gilt auch für die Rate der Oberflächenphononen
    auto t0 = std::chrono::high_resolution_clock::now();
    // must be a double in case of <1 phonons spawns for small cells, small dt or large n_eff
    std::unordered_map<unsigned, std::unordered_map<unsigned, double>> n_phonons;
    unsigned n_spawned=0;

    for(auto it = cell_volumes.begin(); it!= cell_volumes.end() ; ++it) {
        auto [i,j,k] = (*it).first;
        unsigned material_idx = 0;
        for(const auto& material : materials) {
            unsigned branch_idx = 0;
            for(const auto& branch : material.branches) {
                vec3 cellpos = cell_to_pos(i,j,k);
                n_phonons[material_idx][branch_idx] += branch->spatial_spawn_rate(
                        deviational_temp(cellpos), arma::norm(deviational_temp_gradient(cellpos))
                ) / n_eff * (*it).second[material_idx] * dt;
                ++branch_idx;    
            }
            ++material_idx;
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    #ifdef BOLTZI_MEASURE_TIMINGS
    std::cout<<"Spatial first step: "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0)<<std::endl;
    #endif 

    unsigned material_idx = 0;
    for(const auto& material : materials) {
        unsigned branch_idx = 0;
        for(const auto& branch : material.branches) { 
            if(use_turnover) {
                n_phonons[material_idx][branch_idx] += spatial_turnover[material_idx][branch_idx];
                spatial_turnover[material_idx][branch_idx] = util::modf(n_phonons[material_idx][branch_idx]);
            }
            n_spawned+=static_cast<unsigned>(n_phonons[material_idx][branch_idx]);
            auto positions = geometry.random_positions(static_cast<unsigned>(n_phonons[material_idx][branch_idx]), 
                [this, &branch](vec3 pos)->double{
                    return branch->spatial_spawn_rate(deviational_temp(pos), arma::norm(deviational_temp_gradient(pos)));
                });

            for(auto pos : positions) {
                short sign = static_cast<short>((util::random.draw_real() >= .5)*2-1);
                vec3 normal = -arma::normalise(deviational_temp_gradient(pos)) * sign;  // negation here is important
                vec3 e_1 = arma::cross(arma::norm(normal-vec3{0,0,1})>arma::norm(normal-vec3{0,1,0}) ? vec3{0,0,1} : vec3{0,1,0}, normal);
                vec3 e_2 = arma::cross(normal, e_1);
                auto k = branch->draw_k_spatial(deviational_temp(pos), arma::norm(deviational_temp_gradient(pos)), normal, e_1, e_2);
                phonons.add(PhononData{
                    k,
                    pos,
                    material_idx,
                    branch_idx,
                    sign,
                    init_time + util::random.draw_real()*dt,
                    -2
                });
            }
            ++branch_idx;
        }
        ++material_idx;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    #ifdef BOLTZI_MEASURE_TIMINGS
    std::cout<<"Spatial second step: "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1)<<std::endl;
    #endif
    return n_spawned;
}
