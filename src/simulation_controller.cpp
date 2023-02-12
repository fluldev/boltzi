#include "headers/boltzi_types.h"
#include "headers/boundary_current_log.h"
#include <armadillo>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>
#include <atomic>
#include <sstream>
#include <thread>
#include <vector>
#ifdef BOLTZI_MEASURE_TIMINGS
#include <iostream>
#endif

#include "headers/sample.h"
#include "headers/simulation_controller.h"
#include "headers/phonon.h"
#include "headers/util/containers.h"
#include "headers/util/random.h"
#include "headers/util/pool.h"
#ifdef BOLTZI_DASHBOARD
#include "headers/dashboard.h"
#endif
#include "headers/silencer.h"
#include "headers/data_tracer.h"


namespace boltzi{
using config::Silencer;
unsigned SimulationController::advection_step(PhononList::phonon_iterator phonon, double dt, double t, bool log_collisions) {
    if(!phonon->active) return 0;
    double remaining_time = dt;
    unsigned n_coll = 0;
    vec3 vel = working_sample.velocity(*phonon);
    vec3 dir = arma::normalise(vel);
    double move_dist = arma::norm(vel)*remaining_time;

    phonon->last_event_path += move_dist;
    phonon->last_event_time += dt;

    auto coll_event = working_sample.geometry.collision_point(phonon->r, dir, nullptr);
    while(coll_event != nullptr) {
        // TODO: warning in case of two many collisions since this could lead to a infinite loop
        // TODO: only accurate to remaining_time
        double coll_dist = arma::norm(coll_event->position - phonon->r);

        // no remaining collision within timestep
        if(coll_dist>move_dist) break;
        
        remaining_time *= 1-coll_dist/move_dist;
        phonon->r = coll_event->position;

        double reset_means = true;
        //TODO
        // Abhaengigkeit von normale bei kollisionsberechnung hier ignorieren
        // Vorzeichen Strom durch \vec k \cdot \vec n (Oberflaechennormale)
        if(!coll_event->boundary->is_active()) {
            if(log_collisions)
                collision_log.at(coll_event->boundary_idx).transition(
                    t+dt-remaining_time, 
                    (coll_event->negative*2-1)*working_sample.omega(*phonon) * phonon->sign, 
                    phonon->k * phonon->sign,
                    phonon->source_idx
                );  // removed phonons are counted as positive
            reset_means = false;
        }
        else {
            n_coll+=1;

            // check whether phonon gets absorbed
            if(util::random.draw_real()<coll_event->boundary->p_absorb(*phonon, coll_event)) {
                if(log_collisions) 
                    collision_log.at(coll_event->boundary_idx).transition(
                        t+dt-remaining_time, 
                        working_sample.omega(*phonon) * phonon->sign, 
                        phonon->k * phonon->sign,
                        phonon->source_idx
                    );  // removed phonons are counted as positive
                working_sample.phonons.remove(phonon);
                return n_coll;
            }

            if(util::random.draw_real()<coll_event->boundary->p_spec(*phonon, coll_event)) {
                // speculary reflect phonon
                // since this is perfectly elastic it doesn't affect mean free path or time
                reset_means = false;
                phonon->k = coll_event->boundary->reflect(*phonon, coll_event); 
            }
            else {
                // diffusively reflect phonon
                phonon->k = coll_event->boundary->diffusive_draw(*phonon, coll_event);
                #ifdef BOLTZI_COLLISION_POINTS
                collision_tracer.trace(std::make_tuple(
                        coll_event->position[0],
                        coll_event->position[1],
                        coll_event->position[2]
                ));
                #endif
            }
            vel = working_sample.velocity(*phonon);
            dir = arma::normalise(vel);

            if(reset_means) {
                phonon->last_event_path = remaining_time * arma::norm(working_sample.velocity(*phonon));
                phonon->last_event_time = remaining_time;
            }
        }
        move_dist = arma::norm(vel)*std::max(remaining_time, dt*config::SAVE_MOVE_STEP*working_sample.geometry.scale);  // because phonons could get stuck on boundaries otherwise
        coll_event = working_sample.geometry.collision_point(phonon->r, dir, coll_event);
    }
    phonon->r += dir * move_dist;
    return n_coll;
}


bool NonLinearDSMC::collision_step(PhononList::phonon_iterator phonon, double dt, double t) {
    if(!phonon->active) return false;
    const auto [two_phonon_rates, three_phonon_rates] = working_sample.scatter_rates(*phonon, working_sample.get_temperature(phonon->r));
    const double two_phonon_total_tau = std::accumulate(two_phonon_rates.begin(), two_phonon_rates.end(), 0.);
    const double three_phonon_total_tau = std::accumulate(three_phonon_rates.begin(), three_phonon_rates.end(), 0.);
    #ifdef BOLTZI_ALTERNATIVE_SCATTER_PROBABILITY
    if(util::random.draw_real()>dt*(two_phonon_total_tau+three_phonon_total_tau)) return false;  // no scattering
    #else
    if(util::random.draw_real()>(1-std::exp(-dt*(two_phonon_total_tau+three_phonon_total_tau)))) return false;  // no scattering
    #endif
    #ifdef BOLTZI_ALLOW_TRANSITIONS
    auto [_tmp, branch_idx] = working_sample->materials[phonon->material_idx].draw_phonon_branch();
    phonon->branch_idx = branch_idx;
    #endif
    if(util::random.draw_real()*(two_phonon_total_tau+three_phonon_total_tau)<two_phonon_total_tau) { // two phonon process 
        phonon->k = util::random.random_dir() * arma::norm(phonon->k);  
    }
    else {  // three phonon process
        const unsigned idx = util::random.draw_index_weighted(three_phonon_rates);
        double pseudo_temp = working_sample.get_pseudo_temperatures(phonon->r).at(phonon->material_idx).at(phonon->branch_idx)[idx];
        double deviational_temp = working_sample.deviational_temp(phonon->r);
        if(pseudo_temp != deviational_temp)
            phonon->k = working_sample.get_branch(*phonon).three_phonon_processes[idx]->draw_k(
                *phonon,
                working_sample.get_temperature(phonon->r),
                pseudo_temp, 
                deviational_temp 
            );
        else
            // TODO: Unsure if this is the correct way of handling this. 
            phonon->k = working_sample.get_branch(*phonon).three_phonon_processes[idx]->draw_k_linearized(*phonon, deviational_temp);
        phonon->init_time = t;  // technically this is a new phonon
    }
    phonon->last_event_time = dt*util::random.draw_real();
    phonon->last_event_path = phonon->last_event_time * arma::norm(working_sample.velocity(*phonon));
    return true;
}


unsigned NonLinearDSMC::generation_step(double dt, double t) {
    unsigned res = 0;
    if(!no_spatial)
        res+=working_sample.spawn_spatial_phonons(dt, t);
    res+=working_sample.spawn_surface_phonons(dt, t, true, &collision_log);
    return res;
}


void NonLinearDSMC::run(double t1, double t0) {
    for(unsigned i = 0; i<working_sample.geometry.boundaries.size();++i)
        if(working_sample.geometry.boundaries[i]->is_generator() || !working_sample.geometry.boundaries[i]->is_active())
            collision_log.emplace(std::make_pair(i,BoundaryCurrentLog(t1,t0,collision_log_binwidth)));
    working_sample.calc_bin_temps();
    std::vector<std::vector<double>> phonon_graph;
    std::vector<std::vector<double>> temperature_graph;
    volatile std::atomic<unsigned> n_boundary_collisions{0};
    volatile std::atomic<unsigned> n_collisions{0};
    unsigned n_new; 
    auto total_t_0 = std::chrono::high_resolution_clock::now();

    for(double t = t0; t<t1 ; t=t+dt) {
        total_t_0 = std::chrono::high_resolution_clock::now();
        manage_artifacts(t);
        manage_savepoints(t);

        #ifdef BOLTZI_MEASURE_TIMINGS
        std::cout<<"artifact timing: "<<std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-total_t_0)<<std::endl;
        #endif

        n_boundary_collisions = 0;
        n_collisions = 0;
        
        auto timing_0 = std::chrono::high_resolution_clock::now();
        util::pool.add_container_job(
            [this,&n_boundary_collisions,t] (PhononList::phonon_iterator it) {
                n_boundary_collisions+=this->advection_step(it, this->dt, t, true); 
            },
            working_sample.phonons.phonons
        );
        util::pool.exec();
        #ifdef BOLTZI_MEASURE_TIMINGS
        std::cout<<"advection timing: "<<std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-total_t_0)<<std::endl;
        #endif

        // must be parallelized within Sample
        working_sample.calc_bin_temps();

        
        timing_0 = std::chrono::high_resolution_clock::now();
        util::pool.add_container_job(
            [this,&n_collisions,t] (PhononList::phonon_iterator it) {
                n_collisions+=this->collision_step(it, this->dt,t); 
            },
            working_sample.phonons.phonons
        );
        util::pool.exec();
        #ifdef BOLTZI_MEASURE_TIMINGS
        std::cout<<"collision timing: "<<std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-total_t_0)<<std::endl;
        #endif
        
        n_new = generation_step(dt,t);

        #ifdef BOLTZI_MEASURE_TIMINGS
        std::cout<<"generation timing: "<<std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-total_t_0)<<std::endl;
        #endif
        #ifdef BOLTZI_DASHBOARD
        phonon_graph.emplace_back(std::vector<double>{t, static_cast<double>(working_sample.phonons.n_active())});
        temperature_graph.emplace_back(std::vector<double>{t, working_sample.avg_temp()});
        #endif
        if(!Silencer::silent() && !Silencer::less_output()) {
            std::cout<<"\033c";
            #ifdef BOLTZI_DASHBOARD
            dashboard(t, phonon_graph, (t-t0)/dt, (t1-t0)/dt, temperature_graph);
            #endif
        }
        double pts = 0;
        double weight_sum = 0; 
        for(auto& [mat_idx, mat] : working_sample.avg_pseudo_temps())
            for(auto& [branch_idx, branch] : mat) {
                double curweight = working_sample.materials[mat_idx].branches[branch_idx]->branch_weight();
                pts += std::accumulate(branch.begin(),branch.end(),0.0)*curweight / branch.size();
                weight_sum+=curweight;
            }
        pts /= weight_sum;

        Silencer::print(std::stringstream{}
            <<"Step: "<<(t-t0)/dt+1<<"/"<<(t1-t0)/dt<<"\t"
            <<"Phonons: "<<working_sample.phonons.n_active()<<"\t"
            <<"Temperature: "<<working_sample.avg_temp()<<"\t"
            <<"Avg. Pseudo Temperature: "<<pts<<"\t"
            <<"Boundary Collisions: "<<n_boundary_collisions*100. / working_sample.phonons.n_active()<<"%\t"
            <<"Scatter Events: "<<n_collisions*100. / working_sample.phonons.n_active()<<"%\t"
            <<"New Phonons: "<< n_new<<"\n");
        #ifdef BOLTZI_MEASURE_TIMINGS
        std::cout<<"total timing: "<<std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-total_t_0)<<std::endl;
        #endif
    }
    #ifdef BOLTZI_COLLISION_POINTS
    collision_tracer.save();
    #endif
    manage_artifacts(t1);
    manage_savepoints(t1);
    int i = 0;
    for(auto& cl : collision_log) {
        cl.second.save((std::stringstream{}<<boundary_current_log_fname<<"_"<<i<<".boltzi_current").str(), working_sample.n_eff, working_sample.geometry.boundaries[cl.first]->area());
        ++i;
    }
}


void LinearDSMC::run(double t1, double t0) {
    for(unsigned i = 0; i<working_sample.geometry.boundaries.size();++i)
        if(working_sample.geometry.boundaries[i]->is_generator() || !working_sample.geometry.boundaries[i]->is_active())
            collision_log.emplace(std::make_pair(i,BoundaryCurrentLog(t1,t0,collision_log_binwidth)));
    // initialize will be called on demand only
    unsigned size_diff = working_sample.phonons.phonons.size();
    Silencer::print("Spawning spatial phonons.\n");
    working_sample.spawn_spatial_phonons(t1-t0, t0, false);  // spawns them with init_time in between t0 and t1 
    Silencer::print(std::stringstream{}<<"Spawned "<<working_sample.phonons.phonons.size()-size_diff<<" spatial phonons.\n");
    size_diff = working_sample.phonons.phonons.size();
    Silencer::print("Spawning surface phonons.\n");
    //TODO: callback that creates entries in boundary_currents
    //working_sample.spawn_surface_phonons(t1-t0, t0, [this](size_t boundary_idx, double time, double energy){log_boundary_transition(boundary_idx, time, energy);}, false);  // spawns them with init_time in between t0 and t1 
    working_sample.spawn_surface_phonons(t1-t0, t0, false, &collision_log);  // spawns them with init_time in between t0 and t1 
    #ifdef BOLTZI_PHONON_SPAWN_TIMES
    std::ofstream file("phonon_spawn_times.list");
    for(auto& ph : working_sample.phonons.phonons)
        file << ph.init_time << ",";
    file.close();
    #endif
    Silencer::print(std::stringstream{}<<"Spawned "<<working_sample.phonons.phonons.size()-size_diff<<" surface phonons.\n");
    
    // Cloneable template for phonons in pre initialized artifact_phonon_lists.
    // The artifact PhononLists are not allowed to change size due to multithreaded access (size can't be bigger than current amount of phonons).
    PhononData empty_phonon{vec3{0,0,0}, vec3{0,0,0}, 0, 0, 0, 0};
    empty_phonon.active = false;

    std::vector<PhononList> savepoint_phonon_lists(artifacts.size());  // PhononList used by the respective savepoints
    
    for(auto& pl : savepoint_phonon_lists)
        pl.phonons = std::vector<PhononData>(working_sample.phonons.phonons.size(), empty_phonon);

    // PhononList phonons now acts like a job queue.
    util::pool.add_container_job(
        [this,t1,&savepoint_phonon_lists] (PhononList::phonon_iterator phonon){
            double t = phonon->init_time;
            if(std::isnan(t)) // initial phonons
                t = 0;

            #ifdef BOLTZI_FIRST_TRAJECTORY
            std::vector<vec3>* trajectory;
            if(phonon == working_sample.phonons.phonons.begin())
                trajectory = new std::vector<vec3>{};
            #endif

            while(t<t1 && phonon->active) {  // phonon could be deactivated by advection step
                double dt_next_scatter = next_scatter_event(phonon);
                // Positioning before actual scattering is important for correctness.
                #ifdef BOLTZI_FIRST_TRAJECTORY
                if(phonon == working_sample.phonons.phonons.begin()) {
                    trajectory->push_back(phonon->r);
                }
                #endif
                if(dt_next_scatter <= dt_max) {
                    for(unsigned i = 0; i<artifacts.size(); ++i) {
                        if(
                            artifacts[i]->evaluation_time<=t+dt_next_scatter 
                            && artifacts[i]->evaluation_time>=t
                        ) {
                            PhononData bu_phonon = *phonon;
                            // move to correct position in case of big dt_max
                            advection_step(phonon, artifacts[i]->evaluation_time - t, t, false);  
                            artifacts[i]->eval_phonon(*phonon, true);
                            *phonon = bu_phonon;
                        }
                    }
                    for(unsigned i = 0; i<savepoints.size(); ++i) {
                        if(
                            std::get<0>(savepoints[i])<=t+dt_next_scatter 
                            && std::get<0>(savepoints[i])>=t
                        ) {
                            PhononData bu_phonon = *phonon;
                            // move to correct position in case of big dt_max
                            advection_step(phonon, std::get<0>(savepoints[i]) - t, t, false);  
                            savepoint_phonon_lists[i].phonons[phonon - working_sample.phonons.phonons.begin()] = *phonon;
                            *phonon = bu_phonon;
                        }
                    }
                    advection_step(phonon, dt_next_scatter, t, true); 
                    t+=dt_next_scatter;
                    scatter(phonon, t);
                }
                else {
                    for(unsigned i = 0; i<artifacts.size(); ++i) {
                        if(
                            artifacts[i]->evaluation_time<=t+dt_max 
                            && artifacts[i]->evaluation_time>=t
                        ) {
                            PhononData bu_phonon = *phonon;
                            // move to correct position in case of big dt_max
                            advection_step(phonon, artifacts[i]->evaluation_time - t, t, false);  
                            artifacts[i]->eval_phonon(*phonon, true);
                            *phonon = bu_phonon;
                        }
                    }
                    for(unsigned i = 0; i<savepoints.size(); ++i) {
                        if(
                            std::get<0>(savepoints[i])<=t+dt_max 
                            && std::get<0>(savepoints[i])>=t
                        ) {
                            PhononData bu_phonon = *phonon;
                            // move to correct position in case of big dt_max
                            advection_step(phonon, std::get<0>(savepoints[i]) - t, t, false);  
                            savepoint_phonon_lists[i].phonons[phonon - working_sample.phonons.phonons.begin()] = *phonon;
                            *phonon = bu_phonon;
                        }
                    }
                    advection_step(phonon, dt_max, t, true); 
                    t+=dt_max;
                }
            }
            working_sample.phonons.remove(phonon);  // required by informer thread
            #ifdef BOLTZI_FIRST_TRAJECTORY
            if(phonon == working_sample.phonons.phonons.begin()) {
                trajectory->push_back(phonon->r);
                std::ofstream file("trajectory.csv");
                for(auto pos : *trajectory)
                    file<<pos[0]<<","<<pos[1]<<","<<pos[2]<<"\n";
                delete trajectory;
                file.close();
            }
            #endif
        },
        working_sample.phonons.phonons, 1   // only one job per worker
    );
    volatile bool stop_informer = false;
    std::thread informer([&]{
        if(Silencer::silent()) return;
        while(!stop_informer) {
            std::cout<<"\033c";
            std::cout<<"Phonons evaluated: "
            <<working_sample.phonons.phonons.size()-working_sample.phonons.n_active()
            <<" = "<<(working_sample.phonons.phonons.size()-working_sample.phonons.n_active())*100. / working_sample.phonons.phonons.size()<<"%"<<std::endl;
            std::this_thread::sleep_for(std::chrono::milliseconds(config::INFORMER_SLEEP_DURATION));
        }
    });
    
    Silencer::print(std::stringstream{}<<"Starting pool with "<<util::pool.n_jobs()<<" jobs.\n");
    util::pool.exec();
    stop_informer = true;
    informer.join();
    #ifdef BOLTZI_COLLISION_POINTS
    collision_tracer.save();
    #endif
    Silencer::print("Evaluating artifacts.\n");
    for(auto& artifact : artifacts)
        artifact->exec();
    Silencer::print("Evaluating savepoints.\n");
    manage_savepoints(savepoint_phonon_lists);
    int i = 0;
    for(auto& cl : collision_log) {
        cl.second.save((std::stringstream{}<<boundary_current_log_fname<<"_"<<i<<".boltzi_current").str(), working_sample.n_eff, working_sample.geometry.boundaries[cl.first]->area());
        ++i;
    }
}

double LinearDSMC::next_scatter_event(PhononList::phonon_iterator phonon) {
    const auto [two_phonon_rates, three_phonon_rates] = working_sample.scatter_rates(*phonon, working_sample.deviational_temp(phonon->r));
    const double two_phonon_total_tau = std::accumulate(two_phonon_rates.begin(), two_phonon_rates.end(), 0.);
    const double three_phonon_total_tau = std::accumulate(three_phonon_rates.begin(), three_phonon_rates.end(), 0.);
    return -std::log(util::random.draw_real())/(three_phonon_total_tau+two_phonon_total_tau);
}

void LinearDSMC::scatter(PhononList::phonon_iterator phonon, double t) {
    const auto [two_phonon_rates, three_phonon_rates] = working_sample.scatter_rates(*phonon, working_sample.deviational_temp(phonon->r));
    const double two_phonon_total_tau = std::accumulate(two_phonon_rates.begin(), two_phonon_rates.end(), 0.);
    const double three_phonon_total_tau = std::accumulate(three_phonon_rates.begin(), three_phonon_rates.end(), 0.);
    #ifdef BOLTZI_ALLOW_TRANSITIONS
    auto [_tmp, branch_idx] = working_sample->materials[phonon->material_idx].draw_phonon_branch();
    phonon->branch_idx = branch_idx;
    #endif
    if(util::random.draw_real()*(two_phonon_total_tau+three_phonon_total_tau)<two_phonon_total_tau) { // two phonon process 
        phonon->k = util::random.random_dir() * arma::norm(phonon->k);  
    }
    else {  // three phonon process
        const unsigned idx = util::random.draw_index_weighted(three_phonon_rates);
        phonon->k = working_sample.materials[phonon->material_idx].branches[phonon->branch_idx]->three_phonon_processes[idx]->draw_k_linearized(
            *phonon,
            working_sample.deviational_temp(phonon->r)
        );
        phonon->init_time = t;  // technically this is a new phonon
    }
    phonon->last_event_time = 0;
    phonon->last_event_path = 0;
}

} // boltzi
