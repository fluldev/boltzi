#pragma once

#include <cmath>

#include "../material.h"
#include "../util/lindisp.h"


namespace boltzi {
namespace lindisp  {
    class LinearDispBranch : public Branch {
        private:
        const double v_g;
        const double omega_debye;  // omega_debye is diffrent for every branch but the maximal 
                                   // amount of states in reciprocal space is the same for all branches
        
        public:
        LinearDispBranch(const double v_g, const double omega_debye) : v_g(v_g), omega_debye(omega_debye) {}
        LinearDispBranch(
                const std::vector<std::shared_ptr<TwoPhononProcess>> two_phonon_processes, 
                const std::vector<std::shared_ptr<ThreePhononProcess>> three_phonon_processes, 
                const double v_g,
                const double omega_debye
        ) : Branch{two_phonon_processes, three_phonon_processes}, v_g(v_g), omega_debye(omega_debye) {}
        
        double disp_relation(const PhononData& phonon) const override {return v_g * arma::norm(phonon.k);}

        vec3 group_velocity(const PhononData& phonon) const override {return v_g * arma::normalise(phonon.k);}
        
        double dos(double omega) const override {return util::lindisp::dos(omega, v_g);}

        double branch_weight() const override {return util::lindisp::branch_weight(omega_debye, v_g);}  // TODO: check if correct
        
        double avg_velocity() const override {
            return v_g; 
        }
        
        double debye_frequency() const override {
            return omega_debye;
        }

        vec3 draw_k(double temp, double deviational_temp) const override {
            return  util::random.random_dir() 
                    * util::random.accept_reject(
                        util::lindisp::canonical_dist_factory<2>(temp, deviational_temp, v_g), 0, omega_debye,
                        util::lindisp::canonical_dist_max<2>(temp, deviational_temp, v_g)
                    ) / v_g;
        }

        //std::vector<vec3> draw_k_batch(double temp, double deviational_temp, unsigned int n) const override {
        //    std::vector<vec3> res;
        //    util::BinDist dist{
        //        util::lindisp::canonical_dist_factory<2>(temp, deviational_temp, v_g),
        //        0, omega_debye, config::BINDIST_BINS
        //    };
        //    for(int i=0;i<n;++i)
        //        res.push_back(util::random.random_dir() * dist());
        //    return res;
        //}
        
        vec3 draw_k_area(double temp, double deviational_temp, vec3 normal, vec3 x, vec3 y) const override {
            return  util::random.random_dir_area(normal, x, y) 
                    * util::random.accept_reject(
                        util::lindisp::canonical_dist_factory<2>(temp, deviational_temp, v_g), 0, omega_debye,
                        util::lindisp::canonical_dist_max<2>(temp, deviational_temp, v_g)
                    ) / v_g;
        }

        //std::vector<vec3> draw_k_area_batch(double temp, double deviational_temp, vec3 normal, vec3 x, vec3 y, unsigned int n) const override {
        //    std::vector<vec3> res;
        //    util::BinDist dist{
        //        util::lindisp::canonical_dist_factory<2>(temp, deviational_temp, v_g),
        //        0, omega_debye, config::BINDIST_BINS
        //    };
        //    for(int i=0;i<n;++i)
        //        res.push_back(dist() * util::random.random_dir_area(normal, x, y));
        //    return res;
        //}
        
        double omegasum_to_temp(double omegasum) const override {
            omegasum = std::max(0.0, omegasum);
            return util::lindisp::canonical_dist_temp(omegasum, v_g, omega_debye);  
        }
        
         double calc_omegasum(double temp) const override {
            return util::lindisp::canonical_dist_integral<3>(v_g, temp, omega_debye);
        }
        
        double phonon_density(double temp, double deviational_temp) const override {
            return std::abs( 
                util::lindisp::canonical_dist_integral<2>(v_g, temp, omega_debye)
                -util::lindisp::canonical_dist_integral<2>(v_g, deviational_temp, omega_debye)
            );
        }

        double phonon_area_flux(double temp, double deviational_temp) const override {
            return phonon_density(temp, deviational_temp) * v_g / 4;
        }
        
        double spatial_spawn_rate(double deviational_temp, double deviational_temp_gradient) const override {
            return util::lindisp::total_spawnrate_spatial<2>(deviational_temp_gradient, v_g, deviational_temp, omega_debye);
        }
        vec3 draw_k_spatial(double deviational_temp, double deviational_temp_gradient, vec3 normal, vec3 e_1, vec3 e_2) const override {
            return util::random.accept_reject(
                [this, deviational_temp, deviational_temp_gradient] (double omega) -> double { 
                    return util::lindisp::spawnrate_dist_spatial<2>(omega, deviational_temp);
                }, 0, omega_debye,
                util::lindisp::spawnrate_dist_spatial_max<2>(deviational_temp)
            ) * util::random.random_dir_area(normal, e_1, e_2) / v_g;
        }
    };


}
}
