#pragma once

#include <fstream>
#include <functional>
#include <memory> 
#include <mutex>
#include <numeric>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "boltzi_config.h"
#include "phonon.h"
#include "sample.h"
#include "util/containers.h"
#include "util/pool.h"
#include "boundary_current_log.h"
#include "data_tracer.h"

#include "silencer.h"
// Non-linear DSMC:
// Spawn phonons at the correct time in the simulation and evaluate all steps for all phonons for each dt with dt in between two time steps. 
// This implies a call to working_sample->calc_bin_temps() at each dt.
//
// Linear DSMC:
// Populate working_sample->phonons at the beginning.
// Draw evenly distributed spawn time in the simulation timespan for each phonon with dt in between two time steps.
// Evaluate all sim steps in a worker for each phonon independently and in parallel (no call to working_sample->calc_bin_temps()).
// No constant timesteps are needed here, run time until next scatter event or collision.
// Save relevant phonon properties for each timestep.


namespace boltzi {
    using config::Silencer;
    class Phonon {
        private:
        const PhononData& data;
        const Sample& sample;

        public:
        Phonon() = delete;
        Phonon(const PhononData& data, const Sample& sample) : data(data), sample(sample) {}

        double frequency() const {return sample.omega(data);}
        vec3 velocity() const {return sample.velocity(data);}
        vec3 wave_number() const {return data.k;}
        short sign() const {return data.sign;}
        vec3 position() const {return data.r;}
        double init_time() const {return data.init_time;}
        double free_path() const {return data.last_event_path;}
        double free_time() const {return data.last_event_time;}
    };


    class SimulationController {
        protected:
        Sample& working_sample;
        std::vector<std::tuple<double,bool,std::string>> savepoints;
        bool no_spatial;
        const double collision_log_binwidth;
        #ifdef BOLTZI_COLLISION_POINTS
        DataTracer<double, double, double> collision_tracer;
        #endif

        public:
        std::string boundary_current_log_fname{"boundarycurrent"};
        std::unordered_map<unsigned, BoundaryCurrentLog> collision_log;
        struct Artifact {
            bool evaluated = false;
            vec3 cell_unit;
            double cell_volume_unit;
            const std::string artifact_name;
            const double evaluation_time;
            const std::tuple<unsigned, unsigned, unsigned> resolution;
            util::Grid<3, std::unordered_map<unsigned, double>> cell_volumes;
            const Sample* owner;
            virtual void exec() = 0;
            virtual void exec(PhononList& _phonon_list) = 0;
            virtual void eval_phonon(const PhononData& phonon, bool parallel=false) = 0;
            Artifact() = delete;
            Artifact(std::string artifact_name, double evaluation_time, std::tuple<unsigned, unsigned, unsigned> resolution) 
                : artifact_name(artifact_name), evaluation_time(evaluation_time), resolution(resolution), 
                cell_volumes{std::get<0>(resolution), std::get<1>(resolution), std::get<2>(resolution)} {}

            vec3 cell_to_pos(unsigned i, unsigned j, unsigned k) {
                return vec3{
                    owner->geometry.box.first[0] + (i+.5)*cell_unit[0],
                    owner->geometry.box.first[1] + (j+.5)*cell_unit[1],
                    owner->geometry.box.first[2] + (k+.5)*cell_unit[2]
                }; 
            }
            
            std::tuple<unsigned,unsigned,unsigned> pos_to_cell(vec3 r) {
                return std::tuple<unsigned, unsigned, unsigned>(
                    std::get<0>(resolution)*(r[0]-owner->geometry.box.first[0]) / (owner->geometry.box.second[0] - owner->geometry.box.first[0]),
                    std::get<1>(resolution)*(r[1]-owner->geometry.box.first[1]) / (owner->geometry.box.second[1] - owner->geometry.box.first[1]), 
                    std::get<2>(resolution)*(r[2]-owner->geometry.box.first[2]) / (owner->geometry.box.second[2] - owner->geometry.box.first[2])
                ); 
            }


            void init_cells() { 
                cell_unit = (owner->geometry.box.second-owner->geometry.box.first) 
                                / vec3{(double)std::get<0>(resolution), (double)std::get<1>(resolution), (double)std::get<2>(resolution)};
                cell_volume_unit = arma::prod(cell_unit) / config::CELL_OVERLAP_SAMPLES;
                Silencer::print("Initializing artifact cells.\n");
                util::pool.add_container_job(
                    [this] (decltype(cell_volumes)::iterator it) {
                        auto [i,j,k] = (*it).first;
                        for(unsigned n = 0; n<config::CELL_OVERLAP_SAMPLES;++n) {
                            vec3 r = util::random.random_pos_cell(
                                        owner->geometry.box.first + cell_unit 
                                        % vec3{
                                            static_cast<double>(i),
                                            static_cast<double>(j),
                                            static_cast<double>(k)
                                        }, 
                                        cell_unit
                                    );
                            if(owner->geometry.is_inside(r))
                                (*it).second[owner->material_index(r)] += cell_volume_unit;               
                        }
                    },
                cell_volumes);
                util::pool.exec();
            }
        };

        std::vector<std::shared_ptr<Artifact>> artifacts;
        
        template<typename T_res, typename T_part>
        struct ActualArtifact : public Artifact {
            // T_part must support default initialization (to 0).
            // T_res must support overloaded << for filestreams.
            using branch_container = std::unordered_map<unsigned, std::unordered_map<unsigned, T_part>>;
            using phonon_count_container = std::unordered_map<unsigned, std::unordered_map<unsigned, unsigned>>;
            util::Grid<3, branch_container> part_grid{std::get<0>(resolution), std::get<1>(resolution), std::get<2>(resolution)};
            util::Grid<3, phonon_count_container> n_phonons_grid{std::get<0>(resolution), std::get<1>(resolution), std::get<2>(resolution)};

            std::mutex artifact_lock;

            std::function<T_res(branch_container,vec3,phonon_count_container)> consolidation_fn;
            std::function<T_part(const Phonon& phonon)> summing_fn;
            
            const std::string artifact_type;  // must be set by children class


            ActualArtifact(
                std::string artifact_name,
                double t_eval, 
                std::tuple<unsigned, unsigned, unsigned> resolution,
                std::function<T_part(const Phonon&)> summing_fn, 
                std::function<T_res(branch_container, vec3, phonon_count_container)> consolidation_fn,
                std::string artifact_type
            ) : Artifact(artifact_name, t_eval, resolution), summing_fn(summing_fn), consolidation_fn(consolidation_fn), artifact_type(artifact_type){}

            void exec() override;
            void exec(PhononList& phonon_list) override;

            void eval_phonon(const PhononData& phonon, bool parallel = false) override {
                if(!phonon.active) return;
                if(!owner->geometry.is_inside_box(phonon.r)) return;
                auto[i,j,k] = pos_to_cell(phonon.r);
                if(parallel)
                    artifact_lock.lock();
                part_grid.at({i,j,k})[phonon.material_idx][phonon.branch_idx] += summing_fn(Phonon(phonon, *owner));
                n_phonons_grid.at({i,j,k})[phonon.material_idx][phonon.branch_idx] += 1;
                if(parallel)
                    artifact_lock.unlock();
            }

        };

        unsigned advection_step(PhononList::phonon_iterator phonon, double dt, double t, bool log_collisions=true); // independent on type of simulation

        SimulationController() = delete;
        SimulationController(Sample& sample, double collision_log_binwidth, bool no_spatial=false) 
            : working_sample(sample), no_spatial(no_spatial), collision_log_binwidth(collision_log_binwidth) 
            #ifdef BOLTZI_COLLISION_POINTS
            , collision_tracer("collision_points")
            #endif

            {}

        virtual void run(double t, double t0 = 0) = 0;

        void initialize(std::function<double(vec3)> temperature, double t=0) {  // independent on type of simulation
            working_sample.spawn_phonons(temperature, t);
            working_sample.overwrite_temps(temperature);
        }

        void initialize(double temperature, double t) {
            initialize([temperature] (vec3)->double {return temperature;}, t);
        }

        void reset() {
            working_sample.phonons.clear();
        }

        void add_artifact(std::shared_ptr<Artifact> artifact) {
            artifacts.push_back(artifact);
            artifacts.back()->owner = &working_sample;
            artifacts.back()->init_cells();
        }
        
        void add_savepoint(double time, std::string name) {
            savepoints.push_back(std::make_tuple(time,false,name));
        }
    };


    class NonLinearDSMC : public SimulationController {
        const double dt;
        public:
        NonLinearDSMC() = delete;
        NonLinearDSMC(Sample& sample, double dt, double collision_log_binwidth, bool no_spatial=false) 
            : SimulationController(sample, collision_log_binwidth, no_spatial), dt(dt) {}
        void run(double t1, double t0=0) override;
        bool collision_step(PhononList::phonon_iterator phonon, double dt, double t); // dependent on type of simulation
        unsigned generation_step(double dt, double t); // dependent on type of simulation

        void manage_artifacts(double t) {
            for(auto& artifact : artifacts)
                if(!artifact->evaluated && (t>=artifact->evaluation_time)) {
                    artifact->evaluated = true;
                    artifact->exec(working_sample.phonons); 
                }
        }
        void manage_savepoints(double t) {
            for(auto& [t_save, executed, name] : savepoints)
                if(t>=t_save && !executed) {
                    Silencer::print(std::stringstream{}<<"Saving simulation state to file \""<<name<<"\".\n");
                    executed = true;
                    std::stringstream s;
                    s<<t_save;
                    name = name+std::string("_")+s.str();
                    working_sample.phonons.save(name+std::string(".boltzi_state"));
                }
        }
    };
    
    class LinearDSMC : public SimulationController {
        const double dt_max;
        public:
        LinearDSMC() = delete;
        LinearDSMC(Sample& sample, double dt_max, double collision_log_binwidth) 
            : SimulationController(sample, collision_log_binwidth, false), dt_max(dt_max) {}
        void run(double t1, double t0=0) override;
        double next_scatter_event(PhononList::phonon_iterator phonon);
        void scatter(PhononList::phonon_iterator phonon, double t);

        void manage_artifacts(std::vector<PhononList>& artifact_phonon_lists) {
            for(unsigned i = 0; i<artifacts.size(); ++i) {
                artifacts[i]->exec(artifact_phonon_lists[i]);
            }
        }
        void manage_savepoints(std::vector<PhononList>& savepoint_phonon_lists) {
            for(unsigned i = 0; i<savepoints.size(); ++i) {
                std::stringstream s;
                if(std::get<2>(savepoints[i]).empty())
                    s<<std::get<0>(savepoints[i]);
                else
                    s<<std::get<2>(savepoints[i]);
                savepoint_phonon_lists[i].save(s.str()+std::string(".boltzi_state"));
            }
        }
    };

template<typename T_res, typename T_part>
void SimulationController::ActualArtifact<T_res,T_part>::exec() {
    util::Grid<3, T_res> final_grid{std::get<0>(resolution), std::get<1>(resolution), std::get<2>(resolution)};
    util::pool.add_container_job(
        [this] (typename decltype(final_grid)::iterator cell_it) {
            auto [i,j,k] = (*cell_it).first;
            for(auto& [material_idx, mat] : part_grid.at({i,j,k}))
                for(auto& [branch_idx, branch] : mat)
                    branch *= owner->n_eff / cell_volumes.at({i,j,k})[material_idx];
            (*cell_it).second = consolidation_fn(
                part_grid.at({i,j,k}), 
                cell_to_pos(i,j,k), 
                n_phonons_grid.at({i,j,k})
            );
        },
        final_grid
    );
    util::pool.exec();
    std::ostringstream t_str;
    t_str << evaluation_time;
    std::ofstream fdump{artifact_name+"_" + t_str.str() + ".boltzi_artifact"};

    fdump<<artifact_type<<"\n";
    fdump<<evaluation_time<<"\n";
    fdump<<std::get<0>(resolution)<<",";
    fdump<<std::get<1>(resolution)<<",";
    fdump<<std::get<2>(resolution)<<"\n";

    for(auto it = final_grid.begin(); it!=final_grid.end();++it){
        auto [i,j,k] = (*it).first;
        double vol = 0;
        for(auto [material_idx, part_vol] : cell_volumes.at({i,j,k}))
            vol += part_vol;
        if(vol<=config::MIN_ARTIFACT_VOLUME*cell_volume_unit) continue;
        vec3 pos = cell_to_pos(i, j, k) / owner->geometry.scale;
        fdump<<pos[0]<<","<<pos[1]<<","<<pos[2]<<"\n";
        fdump<<vol<<","<<(*it).second<<"\n";   // scale must be reversed here but not for volume
    }
    fdump.close();
}


template<typename T_res, typename T_part>
void SimulationController::ActualArtifact<T_res,T_part>::exec(PhononList& phonon_list) {
    for(auto phonon = phonon_list.phonons.begin();phonon!=phonon_list.phonons.end();++phonon)
        eval_phonon(*phonon);
    exec();
}
}
