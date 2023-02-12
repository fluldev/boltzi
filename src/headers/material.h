#pragma once

#include <memory>
#include <array>
#include <utility>
#include <vector>
#include <random>
#include <cmath>

#include "boltzi_types.h"
#include "boltzi_util.h"
#include "phonon.h"
#include "boltzi_config.h"


namespace boltzi {
    class TwoPhononProcess{
        public:
        virtual double scatter_rate(const PhononData& phonon, double temp) const = 0;
        // TODO:
        // virtual vec3 draw_k_elastic(const PhononData& phonon) const = 0;  // used for two phonon process scatter event (trivial case is just random direction with same magnitude)
    };

    
    class ThreePhononProcess{
        public:
        virtual double scatter_rate(const PhononData& phonon, double temp) const = 0;
        virtual double pseudo_temp(double scatter_sum, double temp) const = 0;
        virtual double scattered_omegasum(double scatter_temp, double temp) const = 0;
        virtual vec3 draw_k(PhononData& phonon, double scatter_temp, double pseudo_temp, double deviational_temp) const = 0;
        virtual vec3 draw_k_linearized(PhononData& phonon, double equ_temp) const = 0;
    };


    class Branch {
        public:
        std::vector<std::shared_ptr<TwoPhononProcess>>two_phonon_processes;
        std::vector<std::shared_ptr<ThreePhononProcess>> three_phonon_processes;
        
        Branch() = default;
        Branch(
                const std::vector<std::shared_ptr<TwoPhononProcess>>& two_phonon_processes, 
                const std::vector<std::shared_ptr<ThreePhononProcess>>& three_phonon_processes
        ) : two_phonon_processes(two_phonon_processes), three_phonon_processes(three_phonon_processes) {}

        virtual double disp_relation(const PhononData& phonon) const = 0;
        virtual vec3 group_velocity(const PhononData& phonon) const = 0;
        virtual double dos(double omega) const = 0;
        virtual double branch_weight() const = 0;
        virtual vec3 draw_k(double temp, double deviational_temp) const = 0;
        //virtual std::vector<vec3> draw_k_batch(double temp, double deviational_temp, unsigned int n) const = 0;
        virtual vec3 draw_k_area(double temp, double deviational_temp, vec3 normal, vec3 x, vec3 y) const = 0;
        //virtual std::vector<vec3> draw_k_area_batch(double temp, double deviational_temp, vec3 normal, vec3 x, vec3 y, unsigned int n) const = 0;
        virtual double omegasum_to_temp(double omegasum) const = 0;
        virtual double calc_omegasum(double temp) const = 0;
        virtual double phonon_density(double temp, double deviational_temp) const = 0;
        virtual double phonon_area_flux(double temp, double deviational_temp) const = 0;
        virtual double avg_velocity() const = 0;
        virtual double spatial_spawn_rate(double deviational_temp, double deviational_temp_gradient) const = 0;
        virtual vec3 draw_k_spatial(double deviational_temp, double deviational_temp_gradient, vec3 normal, vec3 e_1, vec3 e_2) const = 0;
        virtual double debye_frequency() const = 0;
    };


    class Material {
        private:
        using branch_iterator = std::vector<std::shared_ptr<Branch>>::iterator;
        std::vector<double> branch_weights;
        
        public:
        std::vector<std::shared_ptr<Branch>> branches;
        
        Material() = delete;
        Material(const std::vector<std::shared_ptr<Branch>>& branches) : branches(branches) {
            for(auto b : branches) branch_weights.push_back(b->branch_weight());
        }

        const Branch& associated_branch(const PhononData& phonon) const {
            return *(branches.at(phonon.branch_idx)); 
        };

        std::pair<const Branch&, int> draw_phonon_branch() const {
            unsigned d = util::random.draw_index_weighted(branch_weights);
            return std::pair<const Branch&, int>(**(branches.begin() + d), d);
        }
    };
}
