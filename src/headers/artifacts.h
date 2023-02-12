#pragma once

#include "simulation_controller.h"
#include "util/constants.h"
#include "util/lindisp.h"
#include "boltzi_config.h"


namespace boltzi {
namespace artifacts{
    class Temperature : public SimulationController::ActualArtifact<double,double> {
        public:
        Temperature(std::string name, double t_eval, std::tuple<unsigned,unsigned,unsigned> resolution) 
            : ActualArtifact<double, double>(name, t_eval, resolution, 
            [](const Phonon& phonon)->double {return phonon.frequency() * phonon.sign();},
            [this] (ActualArtifact::branch_container s, vec3 cell_pos, ActualArtifact::phonon_count_container n_phonons)->double {
                double temp = 0;
                double total_weight = 0;
                for(auto& [material_idx, mat] : s)
                    for(auto& [branch_idx, branch_sum] : mat) {
                        auto& branch = owner->materials[material_idx].branches[branch_idx];
                        if(n_phonons[material_idx][branch_idx]<config::MINIMUM_TEMPERATURE_PHONONS) {
                            temp += owner->deviational_temp(cell_pos) * config::MINIMUM_TEMPERATURE_PHONONS;
                            total_weight += config::MINIMUM_TEMPERATURE_PHONONS;
                        }
                        else {
                            temp += util::lindisp::canonical_dist_temp(
                                // init volume cells for each material!
                                branch_sum+util::lindisp::canonical_dist_integral<3>(
                                    branch->avg_velocity(), 
                                    owner->deviational_temp(cell_pos), 
                                    branch->debye_frequency()
                                ), 
                                branch->avg_velocity(), branch->debye_frequency()
                            ) * n_phonons[material_idx][branch_idx];
                            total_weight += n_phonons[material_idx][branch_idx];
                        }
                }
                return temp / total_weight;
            }, "temperature") {
            }
    };
    
    class HeatCurrent : public SimulationController::ActualArtifact<vec3,vec3> {
        public:
        HeatCurrent(std::string name, double t_eval, std::tuple<unsigned,unsigned,unsigned> resolution) 
            : ActualArtifact<vec3, vec3>(name, t_eval, resolution, 
            [](const Phonon& phonon)->vec3 {return phonon.frequency() * phonon.sign() * phonon.velocity();},
            [this] (ActualArtifact::branch_container s, vec3 cell_pos, ActualArtifact::phonon_count_container n_phonons)->vec3 {
                vec3 heatcurrent = vec3{0,0,0};
                for(auto& [material_idx, mat] : s)
                    for(auto& [branch_idx, branch_sum] : mat) {
                        auto& branch = owner->materials[material_idx].branches[branch_idx];
                        heatcurrent += branch_sum * util::HBAR;
                    }
                return heatcurrent;
            }, "heatcurrent") {}
    };
    
    class MeanFreePath : public SimulationController::ActualArtifact<double,double> {
        public:
        MeanFreePath(std::string name, double t_eval, std::tuple<unsigned,unsigned,unsigned> resolution) 
            : ActualArtifact<double, double>(name, t_eval, resolution, 
            [](const Phonon& phonon)->double {return phonon.free_path();},
            [this] (ActualArtifact::branch_container s, vec3 cell_pos, ActualArtifact::phonon_count_container n_phonons)->double {
                double meanfreepath = 0;
                for(auto& [material_idx, mat] : s)
                    for(auto& [branch_idx, branch_sum] : mat) {
                        auto [i,j,k] = pos_to_cell(cell_pos);
                        meanfreepath += branch_sum / n_phonons[material_idx][branch_idx] / owner->n_eff * cell_volumes.at({i,j,k})[material_idx];
                    }
                return meanfreepath;
            }, "meanfreepath") {}
    };
    
    class MeanFreeTime : public SimulationController::ActualArtifact<double,double> {
        public:
        MeanFreeTime(std::string name, double t_eval, std::tuple<unsigned,unsigned,unsigned> resolution) 
            : ActualArtifact<double, double>(name, t_eval, resolution, 
            [](const Phonon& phonon)->double {return phonon.free_time();},
            [this] (ActualArtifact::branch_container s, vec3 cell_pos, ActualArtifact::phonon_count_container n_phonons)->double {
                double meanfreetime = 0;
                for(auto& [material_idx, mat] : s)
                    for(auto& [branch_idx, branch_sum] : mat) {
                        auto [i,j,k] = pos_to_cell(cell_pos);
                        meanfreetime += branch_sum / n_phonons[material_idx][branch_idx] / owner->n_eff * cell_volumes.at({i,j,k})[material_idx];
                    }
                return meanfreetime;
            }, "meanfreetime") {}
    };
    class SimulatedPhononCount : public SimulationController::ActualArtifact<double,double> {
        public:
        SimulatedPhononCount(std::string name, double t_eval, std::tuple<unsigned,unsigned,unsigned> resolution) 
            : ActualArtifact<double, double>(name, t_eval, resolution, 
            [](const Phonon& phonon)->double {return 1;},
            [this] (ActualArtifact::branch_container s, vec3 cell_pos, ActualArtifact::phonon_count_container n_phonons)->double {
                double total_n = 0;
                for(auto& [material_idx, mat] : s)
                    for(auto& [branch_idx, branch_sum] : mat) {
                        auto [i,j,k] = pos_to_cell(cell_pos);
                        total_n += branch_sum / owner->n_eff * cell_volumes.at({i,j,k})[material_idx];
                    }
                return total_n;
            }, "simulatedphonons") {}
    };
    class InternalEnergyDensity : public SimulationController::ActualArtifact<double,double> {
        public:
        InternalEnergyDensity(std::string name, double t_eval, std::tuple<unsigned,unsigned,unsigned> resolution) 
            : ActualArtifact<double, double>(name, t_eval, resolution, 
            [](const Phonon& phonon)->double {return phonon.frequency();},
            [this] (ActualArtifact::branch_container s, vec3 cell_pos, ActualArtifact::phonon_count_container n_phonons)->double {
                double total_e = 0;
                for(auto& [material_idx, mat] : s)
                    for(auto& [branch_idx, branch_sum] : mat) {
                        auto [i,j,k] = pos_to_cell(cell_pos);
                        total_e += branch_sum * util::HBAR;
                    }
                return total_e;
            }, "internalenergydensity") {}
    };
}
}
