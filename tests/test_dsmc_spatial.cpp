#include <fstream>
#include <initializer_list>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>

#include <armadillo>

#include "sample_geometry_factory.h"
#include "simulation_controller.h"
#include "scatter_processes/simple_parametric.h"
#include "branches/linear_dispersion_branch.h"
#include "samples/simple_sample.h"
#include "boundaries/specular_vacuum_boundary.h"
#include "boundaries/simple_isothermal_boundary.h"
#include "boltzi_util.h"
#include "artifacts.h"


using namespace boltzi;

int main() {
    double omega_debye = 300e11;
    double n_eff = 2e30;
    double v_g = 1;
    double init_temp = 10;
    double deviational_temp = 11;
    
    SampleGeometry geometry = load_geometry("../trichter.obj");

    std::cout<<"Loaded SampleGeometry."<<std::endl;
    
    std::initializer_list<std::shared_ptr<TwoPhononProcess>> two_scatter_prc = {
        std::shared_ptr<TwoPhononProcess>(new scatter::SimpleParametricTwoPhononProcess<2,0>(15e-26, v_g, omega_debye))
    };
    std::initializer_list<std::shared_ptr<ThreePhononProcess>> three_scatter_prc = {
        std::shared_ptr<ThreePhononProcess>(new scatter::SimpleParametricThreePhononProcess<1,1>(5e-14, v_g, omega_debye))
    };

    std::vector<std::shared_ptr<Branch>> branches{
        std::shared_ptr<Branch>(new lindisp::LinearDispBranch(
            std::vector<std::shared_ptr<TwoPhononProcess>>(two_scatter_prc),
            std::vector<std::shared_ptr<ThreePhononProcess>>(),
            v_g, omega_debye
        ))  
    };

    Material material(branches);

    samples::SimpleSample sample(
            n_eff,  // n_eff
            std::make_tuple(10,1,1),  // n_bins
            geometry,
            material,
            [init_temp] (vec3 r)->double {return init_temp;},
            [] (vec3 r)->vec3 {return vec3{0,0,0};}
    );

    NonLinearDSMC simulation(sample, 5e-3);
    simulation.add_artifact(std::shared_ptr<SimulationController::Artifact>(new artifacts::Temperature(9.9, std::make_tuple(20,10,10), v_g, omega_debye)));
    simulation.add_artifact(std::shared_ptr<SimulationController::Artifact>(new artifacts::Temperature(19.9, std::make_tuple(20,10,10), v_g, omega_debye)));
    simulation.add_artifact(std::shared_ptr<SimulationController::Artifact>(new artifacts::Temperature(29.9, std::make_tuple(20,10,10), v_g, omega_debye)));
    std::cout<<"Registered artifacts."<<std::endl;
    sample.spawn_phonons([init_temp] (vec3)->double {return init_temp;});
    std::cout<<"Spawned initial phonons."<<std::endl;

    simulation.run(10);
}
