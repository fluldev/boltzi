#pragma once

#include <functional>
#include <algorithm>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>
#include <tuple>
#include <functional>
#include <unordered_map>

#include <fstream>

#include "boltzi_exception.h"
#include "boltzi_types.h"
#include "boltzi_util.h"
#include "phonon.h"
#include "material.h"
#include "boundary.h"
#include "boundary_current_log.h"

namespace boltzi {
    class Sample {
        protected:
        const std::tuple<unsigned, unsigned, unsigned> n_bins;  // TODO: should be named n_cells (n_cells is unfortunately used at other places too)
        const vec3 cell_size;
        util::Grid<3, double> cell_temps;
        util::Grid<3, std::unordered_map<unsigned, std::unordered_map<unsigned, double*>>> cell_pseudo_temps;

        std::unordered_map<unsigned, std::unordered_map<unsigned, double>> spatial_turnover;
        std::unordered_map<unsigned, std::unordered_map<unsigned, double>> volume_turnover;
        std::unordered_map<unsigned, double> area_turnover;

        public:
        std::vector<Material> materials;
        const double n_eff;
        util::Grid<3, std::unordered_map<unsigned, double>> cell_volumes;  // volume overlap sample, cell for each material
        std::function<double(vec3 r)> deviational_temp;  // boundaries need this function to calculate propertys at the exact source position
        std::function<vec3(vec3 r)> deviational_temp_gradient;
        SampleGeometry geometry;
        PhononList phonons;

        double avg_temp() const {
            //TODO: Doesn't mean by volume.
            return std::accumulate(cell_temps.get_data().begin(), cell_temps.get_data().end(), 0.) / cell_temps.get_data().size();
        }

        std::unordered_map<unsigned, std::unordered_map<unsigned, std::vector<double>>> avg_pseudo_temps() {
            std::unordered_map<unsigned, std::unordered_map<unsigned, std::vector<double>>> res;
            for(unsigned material_idx = 0; material_idx<materials.size(); ++material_idx)
                for(unsigned branch_idx = 0; branch_idx<materials[material_idx].branches.size(); ++branch_idx)
                    for(unsigned prc_idx = 0; prc_idx<materials[material_idx].branches[branch_idx]->three_phonon_processes.size(); ++prc_idx) {
                        double sum = 0;
                        unsigned n = 0;
                        for(auto it = cell_pseudo_temps.get_data().cbegin(); it!=cell_pseudo_temps.get_data().cend(); ++it)    
                            if(it->contains(material_idx) && it->at(material_idx).contains(branch_idx)) {
                                n+=1;
                                // TODO: Doesn't mean by volume.
                                sum += it->at(material_idx).at(branch_idx)[prc_idx];
                            }
                        res[material_idx][branch_idx].push_back(sum / n);
                    }
            return res;
        }

        void overwrite_temps(std::function<double(vec3)> temperature) {
            for(unsigned i = 0; i<std::get<0>(n_bins); ++i)
                for(unsigned j = 0; j<std::get<1>(n_bins); ++j)
                    for(unsigned k = 0; k<std::get<2>(n_bins); ++k) {
                        vec3 cur_pos = cell_to_pos(i, j, k);
                        cell_temps.at({i,j,k}) = temperature(cur_pos);
                        for(auto& [material_idx, material] : cell_pseudo_temps.at({i,j,k}))
                            for(auto& [branch_idx, branch] : material)
                                for(unsigned prc_idx = 0; prc_idx<materials[material_idx].branches[branch_idx]->three_phonon_processes.size(); ++prc_idx)
                                    branch[prc_idx] = temperature(cur_pos);
                    }
        }

        Sample(
            double n_eff,
            std::tuple<unsigned, unsigned, unsigned> n_bins,
            const SampleGeometry& geometry,
            const std::vector<Material>& materials,
            std::function<double(vec3)> deviational_temp = [] (vec3) -> double {return 0;},
            std::function<vec3(vec3)> deviational_temp_gradient = [] (vec3) -> vec3 {return vec3{0,0,0};}
        );
        virtual ~Sample();

        
        virtual const Material& associated_material(const PhononData& phonon) const = 0;
        virtual const Material& associated_material(vec3 r) const = 0;
        virtual unsigned material_index(vec3 r) const = 0;

        const Branch& get_branch(const PhononData& phonon) const {
            return associated_material(phonon).associated_branch(phonon);
        }

        vec3 velocity(const PhononData& phonon) const { 
            return get_branch(phonon).group_velocity(phonon);
        }

        double omega(const PhononData& phonon) const {
            return get_branch(phonon).disp_relation(phonon);
        }

        double get_temperature(vec3 pos) const {
            // TODO: Interpolate Temperature in between bins and use this function everywhere
            auto [i,j,k] = pos_to_cell(pos);
            return cell_temps.at({i,j,k});
        }

        const std::unordered_map<unsigned,std::unordered_map<unsigned,double*>>& get_pseudo_temperatures(vec3 pos) const {
            auto [i,j,k] = pos_to_cell(pos);
            return cell_pseudo_temps.at({i,j,k});
        }

        std::pair<std::vector<double>,std::vector<double>> scatter_rates(const PhononData& phonon, double temp) const {
            std::vector<double> res1,res2;

            for(auto prc : associated_material(phonon).associated_branch(phonon).two_phonon_processes)
               res1.push_back(prc->scatter_rate(phonon, temp));

            for(auto prc : associated_material(phonon).associated_branch(phonon).three_phonon_processes)
               res2.push_back(prc->scatter_rate(phonon, temp));

            return std::make_pair(res1, res2);
        }
        
        std::tuple<unsigned, unsigned, unsigned> pos_to_cell(vec3 r) const {
            return std::make_tuple(
                static_cast<unsigned int>((r[0]-geometry.box.first[0]) / (geometry.box.second[0]-geometry.box.first[0]) * std::get<0>(n_bins)),
                static_cast<unsigned int>((r[1]-geometry.box.first[1]) / (geometry.box.second[1]-geometry.box.first[1]) * std::get<1>(n_bins)), 
                static_cast<unsigned int>((r[2]-geometry.box.first[2]) / (geometry.box.second[2]-geometry.box.first[2]) * std::get<2>(n_bins)) 
            );
        }

        vec3 cell_to_pos(unsigned i, unsigned j, unsigned k) const {
            return vec3{
                geometry.box.first[0] + (geometry.box.second[0]-geometry.box.first[0]) * (i+.5) / std::get<0>(n_bins),
                geometry.box.first[1] + (geometry.box.second[1]-geometry.box.first[1]) * (j+.5) / std::get<1>(n_bins),
                geometry.box.first[2] + (geometry.box.second[2]-geometry.box.first[2]) * (k+.5) / std::get<2>(n_bins)
            };
        }

        void calc_bin_temps();  // always needed but only for non linear DSMC within the main loop

        unsigned spawn_spatial_phonons(double dt, double init_time, bool use_turnover = true);

        unsigned spawn_phonons(std::function<double(vec3 pos)> init_temp, double init_time, bool use_turnover = true);

        unsigned spawn_surface_phonons(double dt, double init_time, bool use_turnover = true, std::unordered_map<unsigned, BoundaryCurrentLog>* collision_log = nullptr);

        // Must be called by each implementation of Sample AFTER complete construction (last step in constructor)
        // and only for nonlinear dsmc
        void init_cells();
    };
}
